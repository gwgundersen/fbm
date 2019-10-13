/* SRC-MC.C - Markov chain Monte Carlo for source models. */

/* Copyright (c) 1995-2007 by Radford M. Neal 
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

#include "phi.h"
#include "misc.h"
#include "rand.h"
#include "log.h"
#include "mc.h"
#include "src.h"
#include "numin.h"
#include "data.h"
#include "src-data.h"


static int mc_app_set_range (mc_dynamic_state *, int *, int *);


/* LOCAL VARIABLES. */

static int initialize_done = 0;	/* Has this all been set up? */

static src_spec *src;	/* Specification for source model */
static det_spec *det;	/* Specification for detector noise */
static flow_spec *flow;	/* Specification for flow model */

static src_params *params;


/* SET UP REQUIRED RECORD SIZES PRIOR TO GOBBLING RECORDS. */

void mc_app_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{ src_record_sizes(logg);
}


/* INITIALIZE AND SET UP DYNAMIC STATE STRUCTURE. */

void mc_app_initialize
( log_gobbled *logg,	/* Records gobbled up from head and tail of log file */
  mc_dynamic_state *ds	/* Structure holding pointers to dynamical state */
)
{ 
  if (!initialize_done)
  {
    mc_app_set_range_ptr = mc_app_set_range;

    src = logg->data['S'];
    det = logg->data['T'];
    flow = logg->data['F'];

    src_check_specs_present(src,det,flow);

    /* Read data to fit to. */

    data_spec = logg->data['D'];

    if (data_spec && logg->actual_size['D'] !=
                       data_spec_size(data_spec->N_inputs,data_spec->N_targets))
    { fprintf(stderr,"Data specification record is the wrong size!\n");
      exit(1);
    }

    if (data_spec!=0)
    { src_data_read(1,0);
    }

    /* Using existing state if present, or allocating new state if not. */
  
    if (logg->data['q'])
    { if (logg->actual_size['q'] != src_params_length(src->highN))
      { fprintf(stderr,"Bad parameter record (%d %d %d)\n",
                logg->actual_size['q'],src_params_length(src->highN),
                src->highN);
        exit(1);
      }
      params = logg->data['q'];
    }
    else
    { params = chk_alloc (src_params_length(src->highN), sizeof(double));
      src_default_parameters (src,det,flow,params);
    }

    initialize_done = 1;
  }

  ds->dim = src_params_length(src->highN) / sizeof(double);
  ds->q = (double*)params;
  ds->temp_state = 0;

  ds->stepsize = chk_alloc (ds->dim, sizeof(double));
}


/* RESET INITIALIZE_DONE IN PREPARATION FOR NEW LOG FILE. */

void src_mc_cleanup(void)
{
  initialize_done = 0;
}


/* SAVE POSITION AND ANY AUXILIARY PART OF STATE. */

void mc_app_save
( mc_dynamic_state *ds,	/* Current dyanamical state */
  log_file *logf,	/* Log file state structure */
  int index		/* Index of iteration being saved */
)
{
  logf->header.type = 'q';
  logf->header.index = index;
  logf->header.size = src_params_length(src->highN);
  log_file_append(logf,ds->q);
}


/* APPLICATION-SPECIFIC SAMPLING PROCEDURES. */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double a,
  double a2,
  mc_iter *it,
  mc_temp_sched *sch
)
{
  int N, i, j;

  params = (src_params *) ds->q;
  N = (int) params->N0;

  if (strcmp(op,"shuffle")==0)
  { 
    struct src_desc tmp;

    for (i = 0; i<N-1; i++)
    { j = i + rand_int(N-i);
      tmp = params->src[i];
      params->src[i] = params->src[j];
      params->src[j] = tmp;
    }

    return 1;
  }

  else if (strcmp(op,"prior-gen-unused")==0)
  { 
    for (i = N; i<src->highN; i++)
    { src_prior_one_src (params, src, det, flow, i);
    }

    return 1;
  }

  else if (strcmp(op,"shift-intensities")==0)
  { 
    int rep;			/* Number of times to repeat this operation  */
    int ionly;			/* Update intensities only? */

    int N;			/* Number of active sources */
    int i, j;			/* Indexes of the two sources selected */
    double inti, intj;		/* Initial intensities for the two sources */
    double total, propi;	/* Total intensity & proportion for source i */
    double mean[3];		/* Weighted means of coordinates for sources */
    double diff[3];		/* Differences in coordinates of sources */

    double log_prior;		/* Log of prior probability of intensities */

    double slice_point;		/* Variables used for the slice sampling */
    double curr, new;
    double low_bnd, high_bnd;

    int n, c;

    ionly = a<0;
    rep = a<0 ? -a : a>0 ? a : 1;

    N = (int)params->N0;
    if (N<2) 
    { return 1;
    }

    /* Do rep slice sampling updates of intensities for a pair of sources. */

    for (n = rep; n>0; n--)
    {
      it->slice_calls += 1;

      /* Randomly choose two active sources. */

      i = rand_int(N);
      j = rand_int(N-1);
      if (j>=i) j += 1;

      /* Get their intensities, and compute total and proportion of source i. 
         Do nothing if the total intensity is zero. */

      inti = pow(params->src[i].Q,1/src->powQ);
      intj = pow(params->src[j].Q,1/src->powQ);
      total = inti + intj;

      if (total==0) continue;

      propi = inti / total;

      /* Find weighted means and differences of coordinates for the sources. */

      if (!ionly)
      { for (c = 0; c<3; c++)
        { diff[c] = params->src[j].coord[c] - params->src[i].coord[c];
          mean[c] = propi * params->src[i].coord[c] 
                     + (1-propi) * params->src[j].coord[c];
        }
      }

      /* Find energy for current proportions, and also the log of the
         prior density for the intensities in the original (not raised
         to power) form that we're using here. */

      curr = propi;

      if (!ds->know_pot) 
      { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
        ds->know_pot = 1;
        it->slice_evals += 1;
      }

      log_prior = src->powQ==1 ? 0 : (src->powQ-1) * log(propi*(1-propi));

      /* Sample level for slice. */
  
      slice_point = ds->pot_energy - log_prior + rand_exp();

      /* Sample from slice, starting with a full interval for the proportion
         of the total for source i. */

      low_bnd = 0;
      high_bnd = 1;

      for (;;)
      { 
        /* Sample a new proportion, and convert to intensities. */

        new = low_bnd + rand_uniopen() * (high_bnd-low_bnd);
        params->src[i].Q = pow(new*total,src->powQ);
        params->src[j].Q = pow((1-new)*total,src->powQ);

        /* Adjust coordinates to have the same weighted means, keeping 
           differences the same. */

        if (!ionly)
        { for (c = 0; c<3; c++)
          { params->src[i].coord[c] = mean[c] - (1-propi) * diff[c];
            params->src[j].coord[c] = mean[c] + propi * diff[c];
          }
        }

        /* Evaluate energy for sampled point, and also the log of the
           prior density for the intensities in the original (not raised
           to power) form that we're using here. */

        mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
        it->slice_evals += 1;
        ds->know_pot =1;

        log_prior = src->powQ==1 ? 0 : (src->powQ-1) * log(new*(1-new));

        /* Exit if sampled point is within the slice. */
        
        if (ds->pot_energy-log_prior<=slice_point)
        { break;
        }

        /* Otherwise, shrink the interval. */

        if (new<curr)
        { low_bnd = new;
        }
        else
        { high_bnd = new;
        }
      }
    }

    return 1;
  }

  else
  { return 0;
  }
}


/* EVALUATE POTENTIAL ENERGY AND ITS GRADIENT.  Gradient computation isn't
   currently implemented. */

void mc_app_energy
( mc_dynamic_state *ds,	/* Current dyanamical state */
  int N_approx,		/* Number of gradient approximations in use */
  int w_approx,		/* Which approximation to use this time */
  double *energy,	/* Place to store energy, null if not required */
  mc_value *grad	/* Place to store gradient, null if not required */
)
{
  int N, c, i, j;
  double v, e, lw, idf, tf;
  double upper_time;

  tf = ds->temp_state ? ds->temp_state->inv_temp : 1.0;
  params = (src_params *) ds->q;

  N = (int) params->N0;

  if (grad!=0)
  { fprintf(stderr,"Computation of gradient of energy not implemented\n");
    exit(1);
  }

  /* Check for violations of constraints in prior, and account for
     density of stop time given start time. */

  *energy = 0;

  if (flow->type=='t' && (params->U<flow->lowU || params->U>flow->highU)
   || params->log_width<det->log_low_width 
   || params->log_width>det->log_high_width
   || params->inv_df<det->inv_high_df
   || params->inv_df>det->inv_low_df
   || params->N0<src->lowN 
   || params->N0>=src->highN+1
  )
  { *energy = 1e30;
    return;
  }
  
  for (i = 0; i<src->highN; i++)
  { upper_time = params->src[i].start + src->max_duration;
    if (upper_time>src->max_stop)
    { upper_time = src->max_stop;
    }
    if (params->src[i].Q<pow(src->lowQ,src->powQ) 
     || params->src[i].Q>pow(src->highQ,src->powQ) 
     || params->src[i].start<0 || params->src[i].start>src->max_start
     || params->src[i].stop<params->src[i].start 
     || params->src[i].stop>upper_time)
    { *energy = 1e30;
      return;
    }
    else if (src->max_stop<1e30 || src->max_duration<1e30)
    { *energy += log (upper_time - params->src[i].start);
    } 

    for (j = 0; j<3; j++)
    { if (params->src[i].coord[j]<src->low[j] 
       || params->src[i].coord[j]>src->high[j])
      { *energy = 1e30;
        return;
      }
    }
  }

  /* Evaluate energy contribution from fit of predicted values with actual
     measurements. */

  /* In this version, an uncertainty sigma is provided for each
     concentration measurement.  If overall scale of uncertainty is
     unknown, it can be included through noise width.  Otherwise, set
     noise width = 1 (fixed value), ie, log_width = 0. */

  if (tf>0)
  { 
    for (c = 0; c<N_train; c++)
    {
      v = src_total (src, flow, params, train_inputs+4*c);
      e = (v - train_targets[2*c])/train_targets[2*c+1];
      lw = params->log_width;
      idf = params->inv_df;
  
      if (idf==0)
      { *energy -= tf * (log_phi(e/exp(lw)) - lw);
      }
      else
      { *energy -= tf * (- ((1/idf+1)/2) * log(1+e*e*idf/exp(2*lw)) 
                         - lw + lgamma((1/idf+1)/2) - lgamma(0.5/idf) 
                         + sqrt(idf/M_PI));
      }
    }
  }
}


/* SAMPLE FROM DISTRIBUTION AT INVERSE TEMPERATURE OF ZERO.  Returns zero
   if this is not possible. */

int mc_app_zero_gen
( mc_dynamic_state *ds	/* Set to randomly sampled state */
)
{ 
  src_prior_generate ((src_params *) ds, src, det, flow);

  return 1;
}


/* SET STEPSIZES FOR EACH COORDINATE. */

void mc_app_stepsizes
( mc_dynamic_state *ds	/* Set to stepsizes for components of state */
)
{ 
  src_params *p;
  int i, j;

  p = (src_params *) ds->stepsize;

  p->N0 = (double)(src->highN-src->lowN)/4;
  if (p->N0>1) p->N0 = 1;

  for (i = 0; i<src->highN; i++)
  { p->src[i].Q = (pow(src->highQ,src->powQ)-pow(src->lowQ,src->powQ))/10;
    p->src[i].start = src->max_start/10;
    p->src[i].stop = src->max_stop>=1e30 && src->max_duration>=1e30 ? 0 
       : src->max_duration<1e30 ? src->max_duration/10 : src->max_stop/10;
                
    for (j = 0; j<3; j++)
    { p->src[i].coord[j] = (src->high[j]-src->low[j])/10;
    }
  }

  p->log_width = det->log_low_width==det->log_high_width ? 0 : 0.1;
  p->inv_df = det->inv_low_df==det->inv_high_df ? 0 : 0.1;

  p->U = (flow->highU-flow->lowU)/10;
}


/* SET APPLICATION SPECIFIC COORDINATE RANGE. */

static int mc_app_set_range
( mc_dynamic_state *ds,
  int *first,
  int *last
)
{ 
  src_params *params;
  int c;

  if (*first>=-1) abort();

  params = (src_params *) ds->q;

  switch ('A' + (-2-*first))
  { 
    case 'N':  /* Number of sources */
    { *first = 0;
      *last  = 0;
      return 1;
    }

    case 'F':  /* Parameters of flow model */
    { *first = &params->U - &params->N0;
      *last  = &params->U - &params->N0;
      return 1;
    }

    case 'D':  /* Parameters of detector model */
    { *first = &params->log_width - &params->N0;
      *last  = &params->inv_df - &params->N0;
      return 1;
    }

    case 'S':  /* Parameters for sources that exist */
    { *first = &params->src[0].coord[0] - &params->N0;
      *last  = &params->src[(int)params->N0-1].stop -&params->N0;
      return 1;
    }

    default:
    { return 0;
    }
  }

}
