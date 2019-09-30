/* NET-MC.C - Interface between neural network and Markov chain modules. */

/* Copyright (c) 1995 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "misc.h"
#include "rand.h"
#include "log.h"
#include "mc.h"
#include "data.h"
#include "net.h"
#include "net-data.h"


/* NETWORK VARIABLES. */

static initialize_done = 0;	/* Has this all been set up? */

static net_arch *arch;		/* Network architecture */
static net_priors *priors;	/* Network priors */

static net_sigmas sigmas;	/* Hyperparameters for network, auxiliary state
				   for Monte Carlo */

static net_params params;	/* Pointers to parameters, which are position
				   coordinates for dynamical Monte Carlo */

static net_params stepsizes;	/* Pointers to stepsizes */
static net_values seconds;	/* Second derivatives */
static double *train_sumsq;	/* Sums of squared training input values */

static net_values *deriv;	/* Derivatives for training cases */
static net_params grad;		/* Pointers to gradient for network parameters*/


/* PROCEDURES. */

static void gibbs_noise (void);

static void gibbs_unit (net_param *, net_sigma *, net_sigma *, 
                        int, net_prior *);

static void gibbs_conn (net_param *, net_sigma *, net_sigma *, net_sigma *,
                        int, int, net_prior *);

static void gibbs_adjustments (net_sigma *, double, int,
                               net_param *, net_sigma *, double,
                               net_param *, net_sigma *, double, int,
                               int, int *, net_param **, net_sigma **, 
                                 net_prior *, int *);

static double sum_squares (net_param *, net_sigma *, int);

static double cond_sigma (double, double, double, double, int);


/* SET UP REQUIRED RECORD SIZES PRIOR TO GOBBLING RECORDS. */

void mc_app_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{
  logg->req_size['A'] = sizeof (net_arch);
  logg->req_size['P'] = sizeof (net_priors);
  logg->req_size['D'] = sizeof (data_specifications);
}


/* INITIALIZE AND SET UP DYNAMIC STATE STRUCTURE.  Skips some stuff
   if it's already been done, as indicated by the initialize_done
   variable. */

void mc_app_initialize
( log_gobbled *logg,	/* Records gobbled up from head and tail of log file */
  mc_dynamic_state *ds	/* Structure holding pointers to dynamical state */
)
{ 
  net_value *value_block;
  int value_count;
  int i, j;

  if (!initialize_done)
  {
    /* Check that required specification records are present. */
  
    if ((arch = logg->data['A'])==0)
    { fprintf(stderr,"No architecture specification in log file\n");
      exit(1);
    }
  
    if ((priors = logg->data['P'])==0)
    { fprintf(stderr,"No prior specification in log file\n");
      exit(1);
    }
  
    /* Locate existing network, if one exists. */
  
    sigmas.total_sigmas = net_setup_sigma_count(arch);
    params.total_params = net_setup_param_count(arch);
  
    sigmas.sigma_block = logg->data['S'];
    params.param_block = logg->data['W'];
  
    grad.total_params = params.total_params;
  
    if (sigmas.sigma_block!=0 || params.param_block!=0)
    {
      if (sigmas.sigma_block==0 || logg->index['S']!=logg->last_index
       || params.param_block==0 || logg->index['W']!=logg->last_index)
      { fprintf(stderr,
          "Network stored in log file is apparently incomplete\n");
        exit(1);
      }
  
      if (logg->actual_size['S'] != sigmas.total_sigmas*sizeof(net_sigma)
       || logg->actual_size['W'] != params.total_params*sizeof(net_param))
      { fprintf(stderr,"Bad size for network record\n");
        exit(1);
      }
  
      net_setup_sigma_pointers (&sigmas, arch);
      net_setup_param_pointers (&params, arch);
    }
    else
    {
      sigmas.sigma_block = chk_alloc (sigmas.total_sigmas, sizeof (net_sigma));
      params.param_block = chk_alloc (params.total_params, sizeof (net_param));
  
      net_setup_sigma_pointers (&sigmas, arch);
      net_setup_param_pointers (&params, arch);
   
      net_prior_generate (&params, &sigmas, arch, priors, 1, 0);
    }

    /* Set up stepsize structure. */
  
    stepsizes.total_params = params.total_params;
    stepsizes.param_block = chk_alloc (params.total_params, sizeof (net_param));
  
    net_setup_param_pointers (&stepsizes, arch);

    /* Set up value structure. */

    value_count = net_setup_value_count(arch);
    value_block = chk_alloc (value_count, sizeof *value_block);
  
    net_setup_value_pointers (&seconds, value_block, arch);
  
    /* Read training data, if any, and allocate space for derivatives. */
  
    data_spec = logg->data['D'];

    train_sumsq = chk_alloc (arch->N_inputs, sizeof *train_sumsq);
    for (j = 0; j<arch->N_inputs; j++) train_sumsq[j] = 0;
  
    if (data_spec!=0)
    { 
      net_data_read (1, 0, arch);
    
      deriv = chk_alloc (N_train, sizeof *deriv);
    
      value_count = net_setup_value_count(arch);
      value_block = chk_alloc (value_count*N_train, sizeof *value_block);
    
      for (i = 0; i<N_train; i++) 
      { net_setup_value_pointers (&deriv[i], value_block+value_count*i, arch);
      }
    
      for (j = 0; j<arch->N_inputs; j++)
      { for (i = 0; i<N_train; i++)
        { train_sumsq[j] += train_values[i].i[j] * train_values[i].i[j];
        }
      }
    }

    /* Make sure we don't do all this again. */

    initialize_done = 1;
  }

  /* Set up Monte Carlo state structure. */

  ds->aux_dim = sigmas.total_sigmas;
  ds->aux     = sigmas.sigma_block;

  ds->dim = params.total_params;
  ds->q   = params.param_block;
  
  ds->stepsize = stepsizes.param_block;
}


/* RESET INITIALIZE_DONE IN PREPARATION FOR NEW LOG FILE. */

void net_mc_cleanup(void)
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
  logf->header.size = sigmas.total_sigmas * sizeof (net_sigma);
  log_file_append (logf, sigmas.sigma_block);

  logf->header.type = 'W';
  logf->header.index = index;
  logf->header.size = params.total_params * sizeof (net_param);
  log_file_append (logf, params.param_block);
}


/* APPLICATION-SPECIFIC SAMPLING PROCEDURE.  Does gibbs sampling for 
   hyperparameters ("sample-hyper"), or for noise levels ("sample-noise"),
   or for both ("sample-sigmas"). */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double pm,
  mc_iter *it
)
{
  int sample_hyper, sample_noise;
  int l, i;

  sample_hyper = sample_noise = 0;

  if (strcmp(op,"sample-sigmas")==0)
  { sample_hyper = sample_noise = 1;
  }
  else if (strcmp(op,"sample-hyper")==0)
  { sample_hyper = 1;
  }
  else if (strcmp(op,"sample-noise")==0)
  { sample_noise = 1;
  }
  else
  { return 0;
  }

  if (sample_noise && arch->data_model=='R')
  { 
    gibbs_noise();
  }

  if (sample_hyper)
  {
    if (arch->has_ti) gibbs_unit (params.ti, sigmas.ti_cm, 0,
                                  arch->N_inputs, &priors->ti);
  
    for (l = 0; l<arch->N_layers; l++)
    {
      if (l>0)
      { if (arch->has_hh[l-1]) 
        { gibbs_conn (params.hh[l-1], sigmas.hh_cm[l-1], sigmas.hh[l-1], 
                      sigmas.ah[l], arch->N_hidden[l-1], arch->N_hidden[l], 
                      &priors->hh[l-1]);
        }
      }
    
      if (arch->has_ih[l]) 
      { gibbs_conn (params.ih[l], sigmas.ih_cm[l], sigmas.ih[l], 
                    sigmas.ah[l], arch->N_inputs, arch->N_hidden[l], 
                    &priors->ih[l]); 
      }

      if (arch->has_bh[l]) gibbs_unit (params.bh[l], sigmas.bh_cm[l], 
                                       sigmas.ah[l], arch->N_hidden[l], 
                                       &priors->bh[l]);
      
      if (arch->has_th[l]) gibbs_unit (params.th[l], sigmas.th_cm[l], 0,
                                       arch->N_hidden[l], &priors->th[l]);

      if (arch->has_ah[l])
      { gibbs_adjustments (sigmas.ah[l], priors->ah[l], arch->N_hidden[l], 
          arch->has_bh[l] ? params.bh[l] : 0, sigmas.bh_cm[l], 
            priors->bh[l].alpha[1],
          arch->has_ih[l] ? params.ih[l] : 0, sigmas.ih[l], 
            priors->ih[l].alpha[2], arch->N_inputs,
          l>0, &arch->has_hh[l-1], &params.hh[l-1], &sigmas.hh[l-1], 
            &priors->hh[l-1], &arch->N_hidden[l-1]);
      }
  
      if (arch->has_ho[l]) 
      { gibbs_conn (params.ho[l], sigmas.ho_cm[l], sigmas.ho[l], sigmas.ao,
                    arch->N_hidden[l], arch->N_outputs, &priors->ho[l]);
      }
          
    }
  
    if (arch->has_io) 
    { gibbs_conn (params.io, sigmas.io_cm, sigmas.io, sigmas.ao,
                  arch->N_inputs, arch->N_outputs, &priors->io);
    }
  
    if (arch->has_bo) gibbs_unit (params.bo, sigmas.bo_cm, sigmas.ao,
                                  arch->N_outputs, &priors->bo);

    if (arch->has_ao)
    { 
      gibbs_adjustments (sigmas.ao, priors->ao, arch->N_outputs, 
        arch->has_bo ? params.bo : 0, 
          sigmas.bo_cm,
          priors->bo.alpha[1],
        arch->has_io ? params.io : 0,   
          sigmas.io, 
          priors->io.alpha[2], 
          arch->N_inputs,
        arch->N_layers, arch->has_ho, params.ho, sigmas.ho,
          priors->ho, arch->N_hidden);
    }
  }
  
  ds->know_pot  = 0;
  ds->know_grad = 0;

  return 1;
}


/* DO GIBBS SAMPLING FOR NOISE SIGMAS. */

static void gibbs_noise (void)
{
  double nalpha, nprec, sum, d, ps;
  net_prior *pr;
  int i, j;

  for (i = 0; i<N_train; i++) net_func (&train_values[i], 0, arch, &params);

  pr = &priors->noise;

  if (pr->alpha[1]!=0 && pr->alpha[2]==0)
  {
    for (j = 0; j<arch->N_outputs; j++)
    {
      sum = pr->alpha[1] * (*sigmas.noise_cm * *sigmas.noise_cm);
      for (i = 0; i<N_train; i++)
      { d = train_values[i].o[j] - train_targets[i*arch->N_outputs+j];
        sum += d*d;
      }

      nalpha = pr->alpha[1] + N_train;
      nprec = nalpha / sum;

      sigmas.noise[j] = net_pick_sigma (1/sqrt(nprec), nalpha);
    }
  }

  if (pr->alpha[1]!=0 && pr->alpha[2]!=0)
  {
    for (j = 0; j<arch->N_outputs; j++)
    {
      ps = pr->alpha[2] * (sigmas.noise[j] * sigmas.noise[j]);

      sum = 0;
      for (i = 0; i<N_train; i++)
      { d = train_values[i].o[j] - train_targets[i*arch->N_outputs+j];
        sum += rand_gamma((pr->alpha[2]+1)/2) / ((ps+d*d)/2);
      }

      sigmas.noise[j] = cond_sigma (*sigmas.noise_cm, pr->alpha[1],
                                    pr->alpha[2], sum, N_train);
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]==0 && pr->alpha[2]==0)
  {
    sum = pr->alpha[0] * (pr->width * pr->width);
    for (i = 0; i<N_train; i++)
    { for (j = 0; j<arch->N_outputs; j++)
      { d = train_values[i].o[j] - train_targets[i*arch->N_outputs+j];
        sum += d*d;
      }
    }

    nalpha = pr->alpha[0] + N_train * arch->N_outputs;
    nprec = nalpha / sum;
    *sigmas.noise_cm = net_pick_sigma (1/sqrt(nprec), nalpha);

    for (j = 0; j<arch->N_outputs; j++)
    { sigmas.noise[j] = *sigmas.noise_cm;
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]==0 && pr->alpha[2]!=0)
  {
    ps = pr->alpha[2] * (*sigmas.noise_cm * *sigmas.noise_cm);

    sum = 0;
    for (i = 0; i<N_train; i++)
    { for (j = 0; j<arch->N_outputs; j++) 
      { d = train_values[i].o[j] - train_targets[i*arch->N_outputs+j];
        sum += rand_gamma((pr->alpha[2]+1)/2) / ((ps+d*d)/2);
      }
    }

    *sigmas.noise_cm = cond_sigma (pr->width, pr->alpha[0],
                                   pr->alpha[2], sum, arch->N_outputs*N_train);

    for (j = 0; j<arch->N_outputs; j++)
    { sigmas.noise[j] = *sigmas.noise_cm;
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]!=0)
  {
    sum = 0;
    for (j = 0; j<arch->N_outputs; j++) 
    { sum += 1 / (sigmas.noise[j] * sigmas.noise[j]);
    }

    *sigmas.noise_cm = cond_sigma (pr->width, pr->alpha[0],
                                   pr->alpha[1], sum, arch->N_outputs);
  }
}


/* DO GIBBS SAMPLING FOR SIGMA ASSOCIATED WITH GROUP OF UNITS. */

static void gibbs_unit
( net_param *wt,	/* Parameters associated with each unit */
  net_sigma *sg_cm,	/* Common sigma controlling parameter distribution */
  net_sigma *adj,	/* Adjustments for each unit, or zero */
  int n,		/* Number of units */
  net_prior *pr		/* Prior for sigmas */
)
{ 
  double nalpha, nprec, sum, ps, d;
  int i;

  if (pr->alpha[0]!=0 && pr->alpha[1]==0)
  {
    nalpha = pr->alpha[0] + n;

    nprec= nalpha / (pr->alpha[0] * (pr->width*pr->width)
                      + sum_squares(wt,adj,n));

    *sg_cm = net_pick_sigma (1/sqrt(nprec), nalpha);
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]!=0)
  { 
    ps = pr->alpha[1] * (*sg_cm * *sg_cm);

    sum = 0;
    for (i = 0; i<n; i++)
    { d = adj==0 ? wt[i] : wt[i]/adj[i];
      sum += rand_gamma((pr->alpha[1]+1)/2) / ((ps+d*d)/2);
    }

    *sg_cm = cond_sigma (pr->width, pr->alpha[0], pr->alpha[1], sum, n);
  }
  
}


/* DO GIBBS SAMPLING FOR SIGMAS ASSOCIATED WITH GROUP OF CONNECTIONS. */

static void gibbs_conn
( net_param *wt,	/* Weights on connections */
  net_sigma *sg_cm,	/* Common sigma controlling weights */
  net_sigma *sg,	/* Individual sigmas for source units */
  net_sigma *adj,	/* Adjustments for each destination unit, or zero */
  int ns,		/* Number of source units */
  int nd,		/* Number of destination units */
  net_prior *pr		/* Prior for sigmas */
)
{ 
  double width, nalpha, nprec, sum, ps, d;
  int i, j;

  width = net_prior_width_scaled(pr,ns);

  if (pr->alpha[1]!=0 && pr->alpha[2]==0)
  {
    for (i = 0; i<ns; i++)
    {
      nalpha = pr->alpha[1] + nd;
      nprec = nalpha / (pr->alpha[1] * (*sg_cm * *sg_cm)
                         + sum_squares(wt+nd*i,adj,nd));

      sg[i] = net_pick_sigma (1/sqrt(nprec), nalpha);
    }
  }

  if (pr->alpha[1]!=0 && pr->alpha[2]!=0)
  { 
    for (i = 0; i<ns; i++)
    { 
      ps = pr->alpha[2] * (sg[i]*sg[i]);

      sum = 0;        
      for (j = 0; j<nd; j++)
      { d = adj==0 ? wt[nd*i+j] : wt[nd*i+j]/adj[j];
        sum += rand_gamma((pr->alpha[2]+1)/2) / ((ps+d*d)/2);
      }  

      sg[i] = cond_sigma (*sg_cm, pr->alpha[1], pr->alpha[2], sum, nd);
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]==0 && pr->alpha[2]==0)
  {
    nalpha = pr->alpha[0] + ns*nd;

    sum = pr->alpha[0] * (width * width);
    for (i = 0; i<ns; i++)
    { sum += sum_squares(wt+nd*i,adj,nd);
    }
    nprec = nalpha / sum;

    *sg_cm = net_pick_sigma (1/sqrt(nprec), nalpha);

    for (i = 0; i<ns; i++)
    { sg[i] = *sg_cm;
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]==0 && pr->alpha[2]!=0)
  {
    ps = pr->alpha[2] * (*sg_cm * *sg_cm);

    sum = 0;        
    for (i = 0; i<ns; i++)
    { for (j = 0; j<nd; j++)
      { d = adj==0 ? wt[nd*i+j] : wt[nd*i+j]/adj[j];
        sum += rand_gamma((pr->alpha[2]+1)/2) / ((ps+d*d)/2);
      }     
    }

    *sg_cm = cond_sigma (width, pr->alpha[0], pr->alpha[2], sum, ns*nd);

    for (i = 0; i<ns; i++)
    { sg[i] = *sg_cm;
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]!=0)
  {
    sum = 0;
    for (i = 0; i<ns; i++) 
    { sum += 1 / (sg[i] * sg[i]);
    }

    *sg_cm = cond_sigma (width, pr->alpha[0], pr->alpha[1], sum, ns);
  }
}


/* DO GIBBS SAMPLING FOR UNIT ADJUSTMENTS. */

static void gibbs_adjustments
( net_sigma *adj,	/* Adjustments to sample for */
  double alpha,		/* Alpha for adjustments */
  int nd,		/* Number of units with adjustments */
  net_param *b,		/* Biases for destination units, or zero */
  net_sigma *s,		/* Sigma associated with biases */
  double a,		/* Alpha for this sigma */
  net_param *w1,	/* First set of weights, or zero */
  net_sigma *s1,	/* Sigmas associated with first set */
  double a1,		/* Alpha for first set */
  int n1,		/* Number of source units for first set */
  int nrem,		/* Number of remaining weight sets */
  int *has,		/* Whether each remaining set is present */
  net_param **wr,	/* Remaining sets of weights */
  net_sigma **sr,	/* Sigmas associated with remaining sets */
  net_prior *ar,	/* Priors for remaining sets */
  int *nr		/* Numbers of source units in remaining sets */
)
{ 
  double nalpha, nprec, sum, d, ad;
  int r, i, j;

  for (i = 0; i<nd; i++)
  {
    nalpha = alpha;
    sum = alpha;

    ad = adj[i];
  
    if (b!=0)
    { 
      nalpha += 1;

      d = b[i];
      if (a==0)
      { d /= *s;
      }
      else
      { d /= net_pick_sigma (sqrt ((a * (*s * *s) + (d*d)/(ad*ad)) 
                                    / (a+1)), a+1);
      } 
      sum += d*d;
    }
  
    if (w1!=0)
    { 
      nalpha += n1;

      for (j = 0; j<n1; j++)
      { d = w1[nd*j+i];
        if (a1==0)
        { d /= s1[j];
        }
        else
        { d /= net_pick_sigma (sqrt ((a1 * (s1[j] * s1[j]) + (d*d)/(ad*ad)) 
                                       / (a1+1)), a1+1);
        } 
        sum += d*d;
      }
    }

    for (r = 0; r<nrem; r++)
    {
      if (has[r])
      { 
        nalpha += nr[r];

        for (j = 0; j<nr[r]; j++)
        { d = wr[r][nd*j+i];
          if (ar[r].alpha[2]==0)
          { d /= sr[r][j];
          }
          else
          { d /= net_pick_sigma 
                 (sqrt ((ar[r].alpha[2] * (sr[r][j] * sr[r][j]) + (d*d)/(ad*ad))
                          / (ar[r].alpha[2]+1)), ar[r].alpha[2]+1);
          } 
          sum += d*d;
        }
      }
    }
  
    nprec = nalpha / sum;
  
    adj[i] = net_pick_sigma (1/sqrt(nprec), nalpha);
  }
}


/* EVALUATE POTENTIAL ENERGY AND ITS GRADIENT. */

void mc_app_energy
( mc_dynamic_state *ds,	/* Current dynamical state */
  int N_approx,		/* Number of gradient approximations in use */
  int w_approx,		/* Which approximation to use this time */
  double *energy,	/* Place to store energy, null if not required */
  mc_value *gr		/* Place to store gradient, null if not required */
)
{
  double log_prob;
  int i, low, high;

  if (gr && gr!=grad.param_block)
  { grad.param_block = gr;
    net_setup_param_pointers (&grad, arch);
  }

  net_prior_prob (&params, &sigmas, &log_prob, gr ? &grad : 0, arch, priors, 2);

  if (energy) *energy = -log_prob;

  if (-log_prob>=1e30)
  { if (energy) *energy = 1e30;
    if (gr)
    { for (i = 0; i<ds->dim; i++) gr[i] = 0;
    }
    return;
  }

  if (data_spec!=0)
  {
    if (N_approx>1 && gr)
    { for (i = 0; i<ds->dim; i++) gr[i] /= N_approx;
    }

    low  = (N_train * (w_approx-1)) / N_approx;
    high = (N_train * w_approx) / N_approx;

    for (i = (energy ? 0 : low); i < (energy ? N_train : high); i++)
    {
      net_func (&train_values[i], 0, arch, &params);
      net_model_prob (&train_values[i], train_targets + data_spec->N_targets*i,
                      &log_prob, gr ? &deriv[i] : 0, arch, priors, &sigmas, 2);
  
      if (energy) *energy -= log_prob;

      if (gr && i>=low && i<high)
      { net_back (&train_values[i], &deriv[i], arch->has_ti ? -1 : 0,
                  arch, &params);
        net_grad (&grad, &params, &train_values[i], &deriv[i], arch);
      }
    }

    if (N_approx>1 && gr)
    { for (i = 0; i<ds->dim; i++) gr[i] *= N_approx;
    }
  }
}


/* EVALUATE CHANGE IN ENERGY WITH TEMPERING CHANGE.  */

int mc_app_energy_diff
( mc_dynamic_state *ds,	/* Current dyanamical state */
  mc_temp_sched *sch,	/* Tempering schedule */
  int dir,		/* Direction of change */
  double *energy	/* Place to store energy */
)
{
  return 0;
}


/* SET STEPSIZES FOR EACH COORDINATE. */

void mc_app_stepsizes
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
  int i, j, k, l;
  double w;

  /* Compute second derivatives of minus log likelihood for unit values. */

  net_model_max_second (seconds.o, arch, priors, &sigmas);

  for (l = arch->N_layers-1; l>=0; l--)
  { 
    for (i = 0; i<arch->N_hidden[l]; i++)
    { 
      seconds.h[l][i] = 0;

      if (arch->has_ho[l])
      { for (j = 0; j<arch->N_outputs; j++)
        { w = sigmas.ho[l][i];
          if (sigmas.ao!=0) w *= sigmas.ao[j];
          seconds.h[l][i] += (w*w) * seconds.o[j];
        }
      }

      if (l<arch->N_layers-1 && arch->has_hh[l])
      { for (j = 0; j<arch->N_hidden[l+l]; j++)
        { w = sigmas.hh[l][i];
          if (sigmas.ah[l+1]!=0) w *= sigmas.ah[l+1][j];
          seconds.h[l][i] += (w*w) * seconds.s[l+1][j];
        }
      }

      seconds.s[l][i] = seconds.h[l][i];
    }
  }

  if (arch->has_ti)
  { 
    for (i = 0; i<arch->N_inputs; i++)
    { 
      seconds.i[i] = 0;

      if (arch->has_io)
      { for (j = 0; j<arch->N_outputs; j++)
        { w = sigmas.io[i];
          if (sigmas.ao!=0) w *= sigmas.ao[j];
          seconds.i[i] += (w*w) * seconds.o[j];
        }
      }

      for (l = 0; l<arch->N_layers; l++)
      { if (arch->has_ih[l])
        { for (j = 0; j<arch->N_hidden[l]; j++)
          { w = sigmas.ih[l][i];
            if (sigmas.ah[l]!=0) w *= sigmas.ah[l][j];
            seconds.i[i] += (w*w) * seconds.s[l][j];
          }
        }
      }
    }
  }

  /* Initialize stepsize variables to second derivatives of minus log prior. */

  net_prior_max_second (&stepsizes, &sigmas, arch, priors);

  /* Add second derivatives of minus log likelihood to stepsize variables. */

  if (arch->has_ti)
  { for (i = 0; i<arch->N_inputs; i++)
    { stepsizes.ti[i] += N_train * seconds.i[i];
    }
  }

  for (l = 0; l<arch->N_layers; l++)
  {
    if (arch->has_th[l])
    { for (i = 0; i<arch->N_hidden[l]; i++)
      { stepsizes.th[l][i] += N_train * seconds.h[l][i];
      }
    }

    if (arch->has_bh[l])
    { for (j = 0; j<arch->N_hidden[l]; j++)
      { stepsizes.bh[l][j] += N_train * seconds.s[l][j];
      }
    }

    if (arch->has_ih[l])
    { for (i = 0; i<arch->N_inputs; i++)
      { for (j = 0; j<arch->N_hidden[l]; j++)
        { stepsizes.ih [l] [i*arch->N_hidden[l] + j] 
            += train_sumsq[i] * seconds.s[l][j];
        }
      }
    }
 
    if (l<arch->N_layers-1 && arch->has_hh[l])
    { for (i = 0; i<arch->N_hidden[l]; i++)
      { for (j = 0; j<arch->N_hidden[l+1]; j++)
        { stepsizes.hh [l] [i*arch->N_hidden[l+1] + j] 
            += N_train * seconds.s[l+1][j];
        }
      }
    }

    if (arch->has_ho[l])
    { for (i = 0; i<arch->N_hidden[l]; i++)
      { for (j = 0; j<arch->N_outputs; j++)
        { stepsizes.ho [l] [i*arch->N_outputs + j] += N_train * seconds.o[j];
        }
      }
    }
  }

  if (arch->has_io)
  { for (i = 0; i<arch->N_inputs; i++)
    { for (j = 0; j<arch->N_outputs; j++)
      { stepsizes.io [i*arch->N_outputs + j] += train_sumsq[i] * seconds.o[j];
      }
    }
  }

  if (arch->has_bo)
  { for (j = 0; j<arch->N_outputs; j++)
    { stepsizes.bo[j] += N_train * seconds.o[j];
    }
  }

  /* Convert from second derivatives to appropriate stepsizes. */

  for (k = 0; k<ds->dim; k++)
  { ds->stepsize[k] = 1/sqrt(ds->stepsize[k]);
  }
}


/* COMPUTE ADJUSTED SUM OF SQUARES OF A SET OF PARAMETERS. */

static double sum_squares
( net_param *wt,	/* Parameters to compute sum of squares for */
  net_sigma *adj,	/* Adjustments, or zero */
  int n			/* Number of of parameters */
)
{
  double sum_sq, d;
  int i;

  sum_sq = 0;

  if (adj==0)
  { for (i = 0; i<n; i++)
    { d = wt[i];
      sum_sq += d*d;
    }
  }
  else
  { for (i = 0; i<n; i++)
    { d = wt[i] / adj[i];
      sum_sq += d*d;
    }
  }

  return sum_sq;
}


/* SAMPLE FROM CONDITIONAL DISTRIBUTION FOR SIGMA.  Draws a random value 
   from the conditional distribution for a sigma that is defined by
   its top-down prior and by the sum of the lower-level precision values 
   that it controls. */

static double cond_sigma
( double width,		/* Width parameter for top-level prior */
  double alpha0,	/* Alpha for top-level prior */
  double alpha1,	/* Alpha for lower-level prior */
  double sum,		/* Sum of lower-level precisions */
  int n			/* Number of lower-level precision values */
)
{
  static int warned = 0;

  double a, w, v;
  int tries;

  w = 1 / (width * width);
  a = n*alpha1 - alpha0;

  if (a<=0)
  { /* Can't handle this situation */
    abort();
  }

  tries = 0;

  do 
  { 
    v = rand_gamma(a/2) / (sum*alpha1/2);
    tries += 1;

    if (tries>10000 && !warned) 
    { fprintf(stderr,"WARNING: Over 10000 tries in one call of cond_sigma\n");
      warned = 1;
    }

  } while (rand_uniform() > exp(-alpha0/(2*w*v)));

  return sqrt(v);
}
