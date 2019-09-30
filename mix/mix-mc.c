/* MIX-MC.C - Interface between mixture model and Markov chain modules. */

/* Copyright (c) 1995-2004 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided 
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from Radford Neal's
 * web page, and note is made of any changes made to the programs.  The
 * programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "misc.h"
#include "rand.h"
#include "ars.h"
#include "log.h"
#include "mc.h"
#include "data.h"
#include "prior.h"
#include "model.h"
#include "mix.h"
#include "mix-data.h"


/* ADJUSTABLE PARAMETERS. */

#define Max_extras 1000		/* Maximum number of extras components allowed
				   for gibbs-ext-indictors */

/* MIXTURE MODEL VARIABLES. */

static int initialize_done = 0;	/* Has this all been set up? */

static mix_spec *mx;		/* Mixture model specification */
static model_specification *model; /* Data model */

static int N_targets;		/* Number of target values in dataset */

static mix_hypers *hyp;		/* Hyperparameters for mixture model */

static int Max_active;		/* Maximum number of active components */
static int N_active;		/* Number of currently active components */

static int have_indicators;	/* Are component indicators present? */
static int have_offsets;	/* Are component offsets present? */
static int have_noise_SD;	/* Are component noise variables present? */

static short *indicators;	/* Pointer to array of indicator variables */
static double *offsets;		/* Pointer to array of offset parameters */
static double *noise_SD;	/* Pointer to array of noise std. dev. */
static int *freq;		/* Pointer to array of component frequencies */

static double *pr;		/* Array of relative component probabilities */

static int *scount;		/* Count of data points for each component */
static double *smean;		/* Sample means for each component */
static double *svar;		/* Sample variances for each component */


/* PROCEDURES. */

static void gibbs_hypers (void), gibbs_params (void), gibbs_indicators (int),
            gibbs_ext_indicators(int), met_indicators (int, mc_iter *, int);

static double case_prob (int, int);


/* SET UP REQUIRED RECORD SIZES PRIOR TO GOBBLING RECORDS. */

void mc_app_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{ 
  mix_record_sizes(logg);
}


/* INITIALIZE AND SET UP DYNAMIC STATE STRUCTURE.  Skips some stuff
   if it's already been done, as indicated by the initialize_done
   variable.  The real dynamic state structure is currently empty. */

void mc_app_initialize
( log_gobbled *logg,	/* Records gobbled up from head and tail of log file */
  mc_dynamic_state *ds	/* Structure holding pointers to dynamical state */
)
{ 
  int i, t;

  if (!initialize_done)
  {
    /* Check that required specification records are present. */
  
    mx     = logg->data['P'];
    model  = logg->data['M'];

    mix_check_specs_present (mx, 0, model);

    N_targets = mx->N_targets;
  
    /* Locate existing state records. */
  
    hyp = logg->data['S'];

    have_indicators = logg->data['I']!=0;
    have_offsets    = logg->data['O']!=0;
    have_noise_SD = logg->data['N']!=0;
  
    if (hyp==0 && (have_indicators || have_offsets || have_noise_SD))
    { fprintf(stderr,"Missing hyperparameter record\n");
      exit(1);
    }

    if (!have_indicators && (have_offsets || have_noise_SD))
    { fprintf(stderr,"Missing component indicators\n");
      exit(1);
    }

    if (hyp==0) 
    { hyp = chk_alloc (1, sizeof *hyp);
      mix_hyper_init (mx, model, hyp);
    }

    /* Read training data, if any. */
  
    data_spec = logg->data['D'];

    if (data_spec!=0 && model==0)
    { fprintf(stderr,"No model specified for data\n");
      exit(1);
    }

    if (data_spec && logg->actual_size['D'] !=
                       data_spec_size(data_spec->N_inputs,data_spec->N_targets))
    { fprintf(stderr,"Data specification record is the wrong size!\n");
      exit(1);
    }

    if (data_spec!=0) 
    { mix_data_read (1, 0, mx, model);
    }
    else
    { N_train = 0;
    }

    if (have_indicators)
    { if (logg->actual_size['I'] != N_train * sizeof *indicators)
      { fprintf(stderr,"Record of component indicators is wrong size!\n");
        exit(1);
      }
    }

    /* Figure out size of auxiliary stuff. */

    Max_active = mx->N_components==0 ? N_train : mx->N_components;
    N_active = 0;

    if (have_indicators)
    {
      indicators = logg->data['I'];

      for (i = 0; i<N_train; i++)
      { if (indicators[i]>=N_active) 
        { N_active = indicators[i] + 1;
        }
      }

      if (N_active>Max_active)
      { fprintf(stderr,"Garbled record of component indicators!\n");
        exit(1);
      }
 
      if (have_offsets 
       && logg->actual_size['O'] != N_active*N_targets * sizeof *offsets)
      { fprintf(stderr,"Wrong size of record of component offsets!\n");
        exit(1);
      }
 
      if (have_noise_SD 
       && logg->actual_size['N'] != N_active*N_targets * sizeof *noise_SD)
      { fprintf(stderr, 
         "Wrong size of record of component noise standard deviations!\n");
        exit(1);
      }
    }

    /* Allocate space for auxiliary stuff, and copy over existing info. */

    ds->aux_dim = (Max_active+Max_extras) * N_targets * sizeof *offsets;

    if (model!=0 && model->type=='R')
    { ds->aux_dim += (Max_active+Max_extras) * N_targets * sizeof *noise_SD;
    }

    ds->aux = chk_alloc (ds->aux_dim, sizeof (double));

    offsets = ds->aux;

    if (have_offsets)
    { mc_value_copy (offsets, (double*)logg->data['O'], N_active*N_targets);
    }

    if (model!=0 && model->type=='R')
    { noise_SD = offsets + (Max_active+Max_extras) * N_targets;
      if (have_noise_SD)
      { mc_value_copy (noise_SD, (double*)logg->data['N'],
                       N_active*N_targets);
      }
    }

    /* If not already present, and we have some training data, allocate space 
       for indicators and initialize them to indicate the first component. */

    if (N_train>0 && !have_indicators)
    { indicators = chk_alloc (N_train, sizeof *indicators);
      for (i = 0; i<N_train; i++) 
      { indicators[i] = 0;
      }
      N_active = 1;
      have_indicators = 1;
    }

    /* If we don't have parameters, and we have some data, set them all to 
       the corresponding hyperparameters. */

    if (N_train>0 && !have_offsets)
    { for (i = 0; i<N_active; i++)
      { for (t = 0; t<N_targets; t++)
        { offsets[i*N_targets+t] = hyp->mean[t];
        }
      }
      have_offsets = 1;
    }

    if (N_train>0 && model!=0 && model->type=='R' && !have_noise_SD)
    { for (i = 0; i<N_active; i++)
      { for (t = 0; t<N_targets; t++)
        { noise_SD[i*N_targets+t] = hyp->noise[t];
        }
      }
      have_noise_SD = 1;
    }

    /* Allocate space for component frequencies, and initialize them.  Also
       allocate "pr" array. */

    freq = chk_alloc (Max_active+Max_extras, sizeof *freq);
    pr = chk_alloc (Max_active+Max_extras, sizeof *pr);

    mix_freq (indicators, N_train, freq, N_active);

    /* Allocate space for sample counts, means, and variances. */

    scount = chk_alloc (Max_active+1, sizeof *scount);
    smean  = chk_alloc (Max_active+1, sizeof *smean);
    svar   = chk_alloc (Max_active+1, sizeof *svar);

    /* Make sure we don't do all this again. */

    initialize_done = 1;
  }

  /* Set up Monte Carlo state structure.  Everything is null except the
     auxiliary part of the state. */

  ds->dim = 0;
  ds->q   = 0;

  ds->temp_state = 0;
  
  ds->stepsize = 0;
}


/* RESET INITIALIZE_DONE IN PREPARATION FOR NEW LOG FILE. */

void mix_mc_cleanup(void)
{
  initialize_done = 0;
}


/* SAVE POSITION AND AUXILIARY PART OF STATE. */

void mc_app_save
( mc_dynamic_state *ds,	/* Current dyanamical state */
  log_file *logf,	/* Log file state structure */
  int index		/* Index of iteration being saved */
)
{ 
  logf->header.type = 'S';
  logf->header.index = index;
  logf->header.size = sizeof *hyp;
  log_file_append (logf, hyp);

  if (have_indicators)
  { logf->header.type = 'I';
    logf->header.index = index;
    logf->header.size = N_train * sizeof *indicators;
    log_file_append (logf, indicators);
  }

  if (have_offsets)
  { logf->header.type = 'O';
    logf->header.index = index;
    logf->header.size = N_targets*N_active * sizeof *offsets;
    log_file_append (logf, offsets);
  }

  if (have_noise_SD)
  { logf->header.type = 'N';
    logf->header.index = index;
    logf->header.size = N_targets*N_active * sizeof *noise_SD;
    log_file_append (logf, noise_SD);
  }
}


/* APPLICATION-SPECIFIC SAMPLING PROCEDURE. */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double pm,
  double pm2,
  mc_iter *it,
  mc_temp_sched *sch
)
{
  if (strcmp(op,"gibbs-hypers")==0)
  {
    gibbs_hypers();

    return 1;
  }

  else if (strcmp(op,"gibbs-params")==0)
  { 
    if (N_train>0)
    { gibbs_params();
    }

    return 1;
  }

  else if (strcmp(op,"gibbs-indicators")==0)
  { 
    if (mx->N_components==0)
    { fprintf(stderr,
    "The gibbs-indicators operation cannot be used with an infinite number\n");
      fprintf(stderr,
    "  of components.  Use met-indicators instead.\n");
      exit(1);
    }
    
    if (N_train>0)
    { gibbs_indicators(0);
    }

    return 1;
  }

  else if (strcmp(op,"gibbs-ext-indicators")==0)
  { 
    if ((int)pm!=pm || pm<-1)
    { fprintf(stderr,"Invalid parameter for gibbs-ext-indicators\n");
      exit(1);
    }

    if (N_train>0)
    { gibbs_ext_indicators (pm==0 ? 1 : pm);
    }

    return 1;
  }

  else if (strcmp(op,"gibbs1-indicators")==0)
  { 
    if (N_train>0)
    { gibbs_indicators(1);
    }

    return 1;
  }

  else if (strcmp(op,"met-indicators")==0)
  {
    met_indicators (pm==0 ? 1 : pm, it, 0);

    return 1;
  }

  else if (strcmp(op,"met1-indicators")==0)
  {
    met_indicators (pm==0 ? 1 : pm, it, 1);

    return 1;
  }

  return 0;
}


/* EVALUATE POTENTIAL ENERGY AND ITS GRADIENT.  Currently a null procedure,
   since there is no real dynamical state. */

void mc_app_energy
( mc_dynamic_state *ds,	/* Current dynamical state */
  int N_approx,		/* Number of gradient approximations in use */
  int w_approx,		/* Which approximation to use this time */
  double *energy,	/* Place to store energy, null if not required */
  mc_value *gr		/* Place to store gradient, null if not required */
)
{
  if (energy) *energy = 0;

  return;
}


/* SAMPLE FROM DISTRIBUTION AT INVERSE TEMPERATURE OF ZERO.  Returns zero
   if this is not possible. */

int mc_app_zero_gen
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
  return 0;
}


/* SET STEPSIZES FOR EACH COORDINATE.  Nothing to do, since there are
   no coordinates. */

void mc_app_stepsizes
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
}


/* DO GIBBS SAMPLING FOR HYPERPARAMETERS. */

struct gibbs_hypers_info { double w, a, a0, a1, s; };

static double gibbs_hypers_logp (double l, double *d, void *vp)
{ struct gibbs_hypers_info *p = vp;
  double t = exp(l); 
  double v;
  *d = p->a/2 - t*p->a0/(2*p->w) + p->a1*p->s/(2*t);
  v = l*p->a/2 - t*p->a0/(2*p->w) - p->a1*p->s/(2*t);
  return v;
}

static void gibbs_hypers (void)
{ 
  struct gibbs_hypers_info info;
  double a, p, d, m, v, w, r;
  int i, x, t, n;

  /* Sample for hyperparameters associated with each target value. */

  for (t = 0; t<N_targets; t++)
  {
    /* Sample for mean offset for this target. */

    m = 0;
    for (x = 0; x<N_active; x++)
    { m += offsets[x*N_targets+t];
    }

    if (mx->mean_prior.width==0)
    { hyp->mean[t] = 0;
    }
    else
    { p = 1 / (mx->mean_prior.width * mx->mean_prior.width);
      d = 1 / (hyp->SD[t] * hyp->SD[t]);
      hyp->mean[t] = (d*m) / (N_active*d+p)  
                       + sqrt(1/(N_active*d+p)) * rand_gaussian();
    }

    /* Sample for standard deviation of offset for this target. */

    if (mx->SD_prior.alpha[1]!=0)
    { 
      v = 0;
      for (x = 0; x<N_active; x++)
      { v += (offsets[x*N_targets+t] - hyp->mean[t]) 
              * (offsets[x*N_targets+t] - hyp->mean[t]);
      }

      a = mx->SD_prior.alpha[1] + N_active;
      p = a / (mx->SD_prior.alpha[1] * hyp->SD_cm*hyp->SD_cm + v);
      hyp->SD[t] = prior_pick_sigma (1/sqrt(p), a);
    }

    if (model->type=='R')
    {
      /* Sample for noise standard deviation for this target, when the 
         noise standard deviations for all components are linked to it. */

      if (model->noise.alpha[1]!=0 && model->noise.alpha[2]==0)
      {
        v = 0;
        n = 0;
        for (i = 0; i<N_train; i++)
        { x = indicators[i];
          r = train_targets[i*N_targets+t];
          v += (r - offsets[x*N_targets+t]) * (r - offsets[x*N_targets+t]);
          n += 1;
        }

        a = model->noise.alpha[1] + n;
        p = a / (model->noise.alpha[1] * hyp->noise_cm * hyp->noise_cm + v);
      
        hyp->noise[t] = prior_pick_sigma (1/sqrt(p), a);

        for (x = 0; x<N_active; x++)
        { noise_SD[x*N_targets+t] = hyp->noise[t];
        }
      }

      /* Sample for noise standard deviation for this target, based on
         the noise standard deviations for each component. */

      if (model->noise.alpha[1]!=0 && model->noise.alpha[2]!=0)
      { 
        info.w  = 1 / (hyp->noise_cm*hyp->noise_cm);
        info.a  = model->noise.alpha[1] - N_active * model->noise.alpha[2];
        info.a0 = model->noise.alpha[1];
        info.a1 = model->noise.alpha[2];

        info.s  = 0;
        for (x = 0; x<N_active; x++)
        { info.s += 1 / (noise_SD[x*N_targets+t]*noise_SD[x*N_targets+t]);
        }

        hyp->noise[t] = exp(-0.5 * ars(log(info.w), log(1+1/sqrt(info.a0)),
                                       gibbs_hypers_logp, &info));
      }
    }
  }

  /* Sample for common standard deviation for offsets, when per-component
     standard deviations are tied to the common value. */

  if (mx->SD_prior.alpha[0]!=0 && mx->SD_prior.alpha[1]==0)
  { 
    v = 0;
    for (t = 0; t<N_targets; t++)
    { for (x = 0; x<N_active; x++)
      { v += (offsets[x*N_targets+t] - hyp->mean[t]) 
              * (offsets[x*N_targets+t] - hyp->mean[t]);
      }
    }

    w = prior_width_scaled(&mx->SD_prior,mx->N_targets);

    a = mx->SD_prior.alpha[0] + N_active*N_targets;
    p = a / (mx->SD_prior.alpha[0] * w*w + v);
    hyp->SD_cm = prior_pick_sigma (1/sqrt(p), a);

    for (t = 0; t<N_targets; t++) 
    { hyp->SD[t] = hyp->SD_cm;
    }
  }

  /* Sample for common standard deviation for offsets, when per-component
     standard deviations are variable. */

  if (mx->SD_prior.alpha[0]!=0 && mx->SD_prior.alpha[1]!=0)
  {
    w = prior_width_scaled(&mx->SD_prior,mx->N_targets);

    info.w  = 1 / (w*w);
    info.a  = mx->SD_prior.alpha[0] - N_targets * mx->SD_prior.alpha[1];
    info.a0 = mx->SD_prior.alpha[0];
    info.a1 = mx->SD_prior.alpha[1];

    info.s  = 0;
    for (t = 0; t<N_targets; t++)
    { info.s += 1 / (hyp->SD[t]*hyp->SD[t]);
    }

    hyp->SD_cm = exp(-0.5 * ars(log(info.w), log(1+1/sqrt(info.a0)),
                                gibbs_hypers_logp, &info));
  }    

  if (model->type=='R')
  {
    /* Sample for common noise standard deviation, when all lower-level
       noise standard deviations are linked to it. */

    if (model->noise.alpha[0]!=0 
     && model->noise.alpha[1]==0 
     && model->noise.alpha[2]==0)
    {
      w = model->noise.width;

      v = 0;
      n = 0;
      for (t = 0; t<N_targets; t++)
      { for (i = 0; i<N_train; i++)
        { x = indicators[i];
          r = train_targets[i*N_targets+t];
          v += (r - offsets[x*N_targets+t]) * (r - offsets[x*N_targets+t]);
          n += 1;
        }
      }

      a = model->noise.alpha[0] + n;
      p = a / (model->noise.alpha[0] * w*w + v);
      
      hyp->noise_cm = prior_pick_sigma (1/sqrt(p), a);

      for (t = 0; t<N_targets; t++)
      { hyp->noise[t] = hyp->noise_cm;
        for (x = 0; x<N_active; x++)
        { noise_SD[x*N_targets+t] = hyp->noise[t];
        }
      }
    }

    /* Sample for common noise standard deviation, based on the noise
       standard deviations for each target and component. */

    if (model->noise.alpha[0]!=0 
     && model->noise.alpha[1]==0 
     && model->noise.alpha[2]!=0)
    {
      w = model->noise.width;

      info.w  = 1 / (w*w);
      info.a  = model->noise.alpha[0] 
                  - N_targets * N_active * model->noise.alpha[2];
      info.a0 = model->noise.alpha[0];
      info.a1 = model->noise.alpha[2];

      info.s  = 0;
      for (t = 0; t<N_targets; t++)
      { for (x = 0; x<N_active; x++)
        { info.s += 1 / (noise_SD[x*N_targets+t]*noise_SD[x*N_targets+t]);
        }
      }

      hyp->noise_cm = exp(-0.5 * ars(log(info.w), log(1+1/sqrt(info.a0)),
                                     gibbs_hypers_logp, &info));

      for (t = 0; t<N_targets; t++)
      { hyp->noise[t] = hyp->noise_cm;
      }
    }

    /* Sample for common noise standard deviation, based on the noise
       standard deviations for each target. */
  
    if (model->noise.alpha[0]!=0 && model->noise.alpha[1]!=0)
    {
      w = model->noise.width;

      info.w  = 1 / (w*w);
      info.a  = model->noise.alpha[0] - N_targets * model->noise.alpha[1];
      info.a0 = model->noise.alpha[0];
      info.a1 = model->noise.alpha[1];

      info.s  = 0;
      for (t = 0; t<N_targets; t++)
      { info.s += 1 / (hyp->noise[t]*hyp->noise[t]);
      }

      hyp->noise_cm = exp(-0.5 * ars(log(info.w), log(1+1/sqrt(info.a0)),
                                     gibbs_hypers_logp, &info));
    }
  }
}


/* DO GIBBS SAMPLING FOR PARAMETERS OF ACTIVE COMPONENTS.  For real-valued
   data, the offsets and noise standard deviations are updated.  For binary
   data, the offsets (which determine probabilities) are updated, using
   adaptive rejection sampling. */

struct gibbs_params_info { double freq; int n; double pmu, pvar; };

static double gibbs_params_logp (double x, double *deriv, void *in)
{ struct gibbs_params_info *info = in;
  double l, d;
  l = -0.5 * (x-info->pmu)*(x-info->pmu) / info->pvar;
  d = - (x-info->pmu) / info->pvar;
  if (info->freq>0) 
  { l += - info->n * info->freq * log(1+exp(-x));
    d += info->n * info->freq * exp(-x) / (1+exp(-x));
  }
  if (info->freq<1) 
  { l += - info->n * (1-info->freq) * log(1+exp(x));
    d += - info->n * (1-info->freq) * exp(x) / (1+exp(x));
  }
  *deriv = d;
  return l;
}

static void gibbs_params (void)
{ 
  struct gibbs_params_info info;
  double p, d, a, v;
  int x, t, i;

  for (t = 0; t<N_targets; t++)
  {
    /* Calculate sample means for this target. */
  
    for (x = 0; x<N_active; x++)
    { smean[x] = 0;
      scount[x] = 0;
    }
  
    for (i = 0; i<N_train; i++)
    { x = indicators[i];
      v = train_targets[i*N_targets+t];
      smean[x] += v;
      scount[x] += 1;
    }
  
    for (x = 0; x<N_active; x++)
    { smean[x] /= scount[x];
    }
  
    /* Do Gibbs sampling for offsets for this target, for all components. */

    for (x = 0; x<N_active; x++)
    {
      if (model->type=='R')
      {
        p = 1 / (hyp->SD[t]*hyp->SD[t]);
        d = scount[x] / (noise_SD[x*N_targets+t]*noise_SD[x*N_targets+t]);

        offsets[x*N_targets+t] = (p*hyp->mean[t] + d*smean[x]) / (p + d)
                                       + sqrt(1/(p+d)) * rand_gaussian();
      }
      else if (model->type=='B')
      { 
        info.freq = smean[x];
        info.n    = scount[x];
        info.pmu  = hyp->mean[t];
        info.pvar = hyp->SD[t]*hyp->SD[t];

        offsets[x*N_targets+t] = ars (0.0, 1.0, gibbs_params_logp, &info);
      }
    }

    if (model->type=='R' && model->noise.alpha[2]!=0)
    {
      /* Find average squared differences from offsets. */
  
      for (x = 0; x<N_active; x++)
      { svar[x] = 0;
      }
  
      for (i = 0; i<N_train; i++)
      { x = indicators[i];
        v = train_targets[i*N_targets+t];
        svar[x] += (v - offsets[x*N_targets+t]) * (v - offsets[x*N_targets+t]);
      }
  
      for (x = 0; x<N_active; x++)
      { svar[x] /= scount[x];
      }

      /* Do Gibbs sampling for noise_variances, for all components. */

      for (x = 0; x<N_active; x++)
      {
        a = model->noise.alpha[2] + scount[x];
        p = a / (model->noise.alpha[2] * hyp->noise[t] * hyp->noise[t]
                   + scount[x] * svar[x]);
      
        noise_SD[x*N_targets+t] = prior_pick_sigma (1/sqrt(p), a);
      }
    }
  }
}


/* DO GIBBS SAMPLING FOR COMPONENT INDICATORS.  The parameter indicates
   whether this is a "gibbs-indicators" or "gibbs1-indicators" operation.
   The "gibbs-indicators" operation is possible only if the number of 
   components in the model is finite. */

static void gibbs_indicators
( int gibbs1
)
{
  double c, m;
  int i, x, t;

  if (!gibbs1 && mx->N_components==0) abort();

  /* For "gibbs-indicators", extend number of active components to the 
     maximum, sampling the parameters of the previously unused ones from 
     the prior. */

  if (!gibbs1)
  {
    while (N_active < mx->N_components)
    {
      freq[N_active] = 0;
  
      for (t = 0; t<N_targets; t++)
      {
        offsets[N_active*N_targets+t] = 
          hyp->mean[t] + hyp->SD[t]*rand_gaussian();
  
        if (model->type=='R')
        { noise_SD[N_active*N_targets+t] = 
            prior_pick_sigma (hyp->noise[t], model->noise.alpha[2]);
        }
      }
  
      N_active += 1;
    }  
  }

  /* Do Gibbs sampling for indicators. */

  if (mx->N_components==0)
  { c = 0;
  }
  else
  { c = hyp->con;
    if (mx->con_prior.scale) c /= mx->N_components;
  }

  for (i = 0; i<N_train; i++)
  {
    if (gibbs1 && freq[indicators[i]]==1) 
    { continue;
    }

    freq[indicators[i]] -= 1;

    for (x = 0; x<N_active; x++)
    { pr[x] = log(freq[x]+c) + case_prob(x,i);
    }

    m = pr[0];
    for (x = 1; x<N_active; x++)
    { if (pr[x]>m) 
      { m = pr[x]; 
      }
    }

    for (x = 0; x<N_active; x++)
    { pr[x] = exp(pr[x]-m);
    }

    indicators[i] = rand_pickd(pr,N_active);

    freq[indicators[i]] += 1;
  }
 
  /* Sort components in decreasing order by frequency, and get rid of 
     inactive ones. */

  mix_sort (indicators, N_train, freq, N_active, offsets, noise_SD, N_targets);

  while (freq[N_active-1]==0)
  { if (N_active==1) abort();
    N_active -= 1;
  }
}


/* DO GIBBS SAMPLING FOR INDICATORS USING EXTRA COMPONENTS.  The parameter 
   is the number of extra components, or -1 for the "no gaps" method. */

static void gibbs_ext_indicators
( int nx		/* Number of extra components */
)
{
  int i, x, n, t;
  double c, r, m;
  int no_gaps;

  if (nx>Max_extras)
  { fprintf (stderr, 
      "Extra components for gibbs-ext-indicators exceeds maximum (%d)\n",
      Max_extras);
    exit(1);
  }

  if (nx==-1)
  { nx = 1;
    no_gaps = 1;
  }
  else
  { no_gaps = 0;
  }

  /* Do Gibbs sampling for indicators. */

  if (mx->N_components==0)
  { c = 0;
  }
  else
  { c = hyp->con;
    if (mx->con_prior.scale) c /= mx->N_components;
  }

  for (i = 0; i<N_train; i++)
  {
    /* Do nothing when the "no gaps" method would. */

    if (no_gaps && freq[indicators[i]]==1 
     && rand_uniform()<(N_active-1.0)/N_active) 
    { continue;
    }

    freq[indicators[i]] -= 1;

    /* Draw components from the prior as required. */

    n = freq[indicators[i]]==0 ? N_active+nx-1 
         : N_active==mx->N_components ? N_active
         : N_active+nx;

    for (x = N_active; x<n; x++)
    {
      freq[x] = 0;
  
      for (t = 0; t<N_targets; t++)
      {
        offsets[x*N_targets+t] = 
          hyp->mean[t] + hyp->SD[t]*rand_gaussian();
  
        if (model->type=='R')
        { noise_SD[x*N_targets+t] = 
            prior_pick_sigma (hyp->noise[t], model->noise.alpha[2]);
        }
      }
    }  

    /* Pick a component. */

    if (mx->N_components==0)
    { r = hyp->con / nx;
    }
    else
    { r = c * (mx->N_components - N_active + (freq[indicators[i]]==0)) / nx;
    }

    if (no_gaps) r /= (N_active - (freq[indicators[i]]==0) + 1);

    for (x = 0; x<n; x++)
    { pr[x] = (freq[x]==0 ? log(r) : log(freq[x]+c)) + case_prob(x,i);
    }

    m = pr[0];
    for (x = 1; x<n; x++)
    { if (pr[x]>m) 
      { m = pr[x]; 
      }
    }

    for (x = 0; x<n; x++)
    { pr[x] = exp(pr[x]-m);
    }

    x = rand_pickd(pr,n);

    /* Use the component picked. */

    if (x<N_active)
    { indicators[i] = x;
    }
    else
    { if (x!=N_active)
      { for (t = 0; t<N_targets; t++)
        { offsets[N_active*N_targets+t] = offsets[x*N_targets+t];
          if (model->type=='R')
          { noise_SD[N_active*N_targets+t] = noise_SD[x*N_targets+t];
          }
        }
      }
      indicators[i] = N_active;
      N_active += 1;
    }

    freq[indicators[i]] += 1;

    /* Sort by frequency and eliminate unused components. */

    mix_sort(indicators, N_train, freq, N_active, offsets, noise_SD, N_targets);

    while (freq[N_active-1]==0)
    { if (N_active==1) abort();
      N_active -= 1;
    }
  }
}


/* DO METROPOLIS UPDATES FOR COMPONENT INDICATORS. */

static void met_indicators
( int n_updates,	/* Number of updates to do for each indicator */
  mc_iter *it,		/* Iteration record */
  int met1		/* Is this a "met1-indicators" operation? */
)
{ 
  int Di, singleton, E, N;
  double c, lp, np, dp, U;
  int i, x, t, u;

  N = mx->N_components;

  c = hyp->con;
  if (N!=0 && mx->con_prior.scale) c /= N;

  /* Do Metropolis-Hastings updates for each indicator. */

  Di = rand_int(N_train); /* Selects a training case to record delta for */

  for (i = 0; i<N_train; i++)
  {
    lp = case_prob (indicators[i], i);

    for (u = 0; u<n_updates; u++)
    {
      E = N - N_active;

      /* Pick from the distribution for this indicator based on frequencies and
         the concentration parameter.  For a "met-indicators" operation, all 
         components are possible, including a new one.  For "met1-indicators",
         the possible components are determined by whether it's a singleton. */

      if (met1)
      { 
        singleton = freq[indicators[i]]==1;

        if (singleton)
        { 
          for (x = 0; x<N_active; x++)
          { if (x==indicators[i] || freq[x]==0)
            { pr[x] = 0;
            }
            else
            { pr[x] = freq[x];
              if (N!=0) pr[x] += c;
            }
          }
  
          x = rand_pickd(pr,N_active);
        }
        else
        { 
          x = N!=0 && E==0 ? indicators[i] : N_active;
        }
      }
      else  
      { 
        for (x = 0; x<N_active; x++)
        { pr[x] = freq[x];
          if (N!=0) pr[x] += c;
        }
  
        pr[indicators[i]] -= 1;
  
        pr[N_active] = N==0 ? c : c*E;
  
        x = rand_pickd(pr,N_active+1);
      }
  
      /* Generate parameters from the prior if we picked a new component. */
  
      if (x==N_active)
      {
        for (t = 0; t<N_targets; t++)
        {
          offsets[N_active*N_targets+t] 
            = hyp->mean[t] + hyp->SD[t]*rand_gaussian();
  
          if (model->type=='R')
          { noise_SD[N_active*N_targets+t] = 
              prior_pick_sigma (hyp->noise[t], model->noise.alpha[2]);
          }
        }
      }  
  
      /* Decide whether to accept the switch to new component, and if so, update
         things as required. */
  
      np = x==indicators[i] ? lp : case_prob (x, i);
      dp = lp - np;

      if (met1)
      { if (N==0)
        { dp += singleton ? log(c/(N_train-1)) : -log(c/(N_train-1));
        }
        else if (singleton || E>0)
        { dp += singleton ?  log(c*(E+1)/(N_train-1+c*(N_active-1))) 
                          : -log(c*E/(N_train-1+c*N_active));
        }
      }

      it->proposals += 1;
      if (i==Di) it->delta = dp;

      U = rand_uniform();

      if (dp<0 || U<exp(-dp))
      { 
        if (x==N_active)
        { N_active += 1;
          if (N!=0 && N_active>N) abort();
          freq[x] = 0;
        }
  
        freq[indicators[i]] -= 1;
        indicators[i] = x;
        freq[indicators[i]] += 1;
  
        lp = np;

        mix_sort (indicators, N_train, freq, N_active, offsets, 
                  noise_SD, N_targets);
  
        while (freq[N_active-1]==0)
        { if (N_active==1) abort();
          N_active -= 1;
        }

        if (i==Di) it->move_point = 1;
      }
      else
      { if (i==Di) it->move_point = 0;
        it->rejects += 1;
      }
    }
  }
}


/* COMPUTE LOG PROBABILITY OF CASE FOR A COMPONENT.  Computes the log
   of the probability of a particular training case, absent a constant
   term, with respect to the distribution defined by a particular
   component.  The computation accesses various global variables. */

static double case_prob
( int x,		/* Index of component */
  int i			/* Index of training case */
)
{ 
  double lp, v, n, o;
  int t;
  
  lp = 0;

  for (t = 0; t<N_targets; t++)
  {
    v = train_targets[i*N_targets+t];
    o = offsets[x*N_targets+t];

    if (model->type=='R')
    { n = noise_SD[x*N_targets+t];
      lp -= log(n);
      lp -= (v-o)*(v-o)/(2*n*n);
    }
   
    else if (model->type=='B')
    { lp -= v==0 ? log(1+exp(o)) : log(1+exp(-o));
    }
  }

  return lp;
}
