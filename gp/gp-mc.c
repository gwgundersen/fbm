/* GP-MC.C - Interface between Gaussian process and Markov chain modules. */

/* Copyright (c) 1996 by Radford M. Neal 
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
#include "matrix.h"
#include "rand.h"
#include "ars.h"
#include "log.h"
#include "mc.h"
#include "data.h"
#include "prior.h"
#include "model.h"
#include "gp.h"
#include "gp-data.h"


/* ADJUSTABLE PARAMETER CONTROLLING OPTIONAL MEMORY USAGE.  This adjustable
   parameter sets the maximum amount of memory that will be used to store
   intermediate results (specifically, terms in the covariances), in order to 
   avoid recomputing them.  Its setting has no effect on the results, but 
   setting it higher may reduce computation time, at the cost of memory usage.*/

#define Max_optional_memory 10000000  /* Max. optional memory usage, in bytes */


/* GAUSSIAN PROCESS VARIABLES. */

static initialize_done = 0;	/* Has this all been set up? */

static gp_spec *gp;		/* Gaussian process specification */
static model_specification *model; /* Data model */
static model_survival *surv;	/* Hazard type for survival model */

static gp_hypers hypers;	/* Hyperparameters for Gaussian process, primary
				   state for Monte Carlo */

static gp_hypers stepsizes;	/* Pointers to stepsizes */
static gp_hypers grad;		/* Pointers to derivatives */

static int have_values;		/* Are case-by-case latent values present? */
static int have_variances;	/* Are case-by-case noise variances present? */

static double *aux;		/* Space for auxiliary part of state (latent
				   values and noise variances for each case) */

static double *latent_values;	/* Pointer into aux for latent values */
static double *noise_variances;	/* Pointer into aux for noise variances */

static double *train_cov;	/* Covariance matrix for training points */
static double *diff_cov;	/* Derivative of covariance matrix for trn pts*/

static double *exp_cov[Max_exp_parts]; /* Exponential terms in covariance,
			          saved to speed up computation of diff_cov */

static double *latent_cov;	/* Covariance for latent values at trn points
				    - same storage as exp_cov[0]. */

static double *reg_matrix;	/* Matrix of regression coefficients - same
				   storage as exp_cov[1]. */

static double *post_cov;	/* Posterior cov. for latent values underlying
				   training points - same storage as diff_cov */

static double *post_mean;	/* Posterior mean of latent values at trn pts */

static double *scr1, *scr2;	/* Scratch vectors */


/* PROCEDURES. */

static void scan_values(int), sample_values(void), sample_variances(void);


/* SET UP REQUIRED RECORD SIZES PRIOR TO GOBBLING RECORDS. */

void mc_app_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{ 
  gp_record_sizes(logg);
}


/* INITIALIZE AND SET UP DYNAMIC STATE STRUCTURE.  Skips some stuff
   if it's already been done, as indicated by the initialize_done
   variable. */

void mc_app_initialize
( log_gobbled *logg,	/* Records gobbled up from head and tail of log file */
  mc_dynamic_state *ds	/* Structure holding pointers to dynamical state */
)
{ 
  int optional_left;
  int i, j, l;

  if (!initialize_done)
  {
    /* Check that required specification records are present. */
  
    gp     = logg->data['P'];
    model  = logg->data['M'];
    surv   = logg->data['V'];

    gp_check_specs_present(gp,0,model,surv);
  
    /* Locate existing state records, if they exist. */
  
    hypers.total_hypers = gp_hyper_count(gp,model);
    hypers.hyper_block = logg->data['S'];

    have_values = logg->data['F']!=0;
    have_variances = logg->data['N']!=0;
  
    if (hypers.hyper_block!=0 || have_values || have_variances)
    {
      if (hypers.hyper_block==0)
      { fprintf(stderr,"Missing hyperparameter record\n");
        exit(1);
      }

      if (logg->actual_size['S'] != hypers.total_hypers*sizeof(double))
      { fprintf(stderr,"Bad size for hyperparameter record\n");
        exit(1);
      }

      gp_hyper_pointers (&hypers, gp, model);
    }
    else
    {
      hypers.hyper_block = chk_alloc (hypers.total_hypers, sizeof (double));
      gp_hyper_pointers (&hypers, gp, model);
   
      gp_prior_generate (&hypers, gp, model, 1, 0, 0);
    }

    /* Set up stepsize structure. */
  
    stepsizes.total_hypers = hypers.total_hypers;
    stepsizes.hyper_block = chk_alloc (hypers.total_hypers, sizeof (double));
  
    gp_hyper_pointers (&stepsizes, gp, model);

    /* Read training data, if any, and set up auxiliary state. */
  
    data_spec = logg->data['D'];

    if (data_spec!=0) 
    { gp_data_read (1, 0, gp, model, surv);
    }
    else
    { N_train = 0;
    }

    aux = chk_alloc (2*gp->N_outputs*N_train, sizeof (double));

    latent_values = aux;
    noise_variances = aux + gp->N_outputs*N_train;

    if (have_values)
    { 
      if (logg->actual_size['F']!=gp->N_outputs*N_train*sizeof(double))
      { fprintf (stderr,
          "Bad size for latent values record (is %d, should be %d)\n",
          logg->actual_size['F'], gp->N_outputs*N_train*sizeof(double));
        exit(1);
      }

      for (i = 0; i<N_train; i++)
      { for (j = 0; j<gp->N_outputs; j++)
        { latent_values [gp->N_outputs*i + j] 
            = ((double*) logg->data['F']) [gp->N_outputs*i + j];
        }
      }
    }

    if (have_variances)
    { 
      if (logg->actual_size['N']!=gp->N_outputs*N_train*sizeof(double))
      { fprintf (stderr,
          "Bad size for noise variances record (is %d, should be %d)\n",
          logg->actual_size['N'], gp->N_outputs*N_train*sizeof(double));
        exit(1);
      }

      for (i = 0; i<N_train; i++)
      { for (j = 0; j<gp->N_outputs; j++)
        { noise_variances [gp->N_outputs*i + j] 
            = ((double*) logg->data['N']) [gp->N_outputs*i + j];
        }
      }
    }

    if (model!=0 && model->type=='R' && model->noise.alpha[2]!=0
     && !have_variances)
    { for (i = 0; i<N_train; i++)
      { for (j = 0; j<gp->N_outputs; j++)
        { noise_variances [gp->N_outputs*i + j] = 1.0;
        }
      }
      have_variances = 1;
    }

    /* Allocate space for covariance and derivative matrices, and scratch. */

    train_cov  = chk_alloc (N_train*N_train, sizeof(double));
    diff_cov   = chk_alloc (N_train*N_train, sizeof(double));
    latent_cov = chk_alloc (N_train*N_train, sizeof(double));
    reg_matrix = chk_alloc (N_train*N_train, sizeof(double));

    post_cov  = diff_cov;

    optional_left = Max_optional_memory;

    for (l = 0; l<gp->N_exp_parts; l++)
    { 
      if (l==0)
      { exp_cov[l] = latent_cov;
      }
      else if (l==1)
      { exp_cov[l] = reg_matrix;
      }
      else if (N_train*N_train*sizeof(double)<=optional_left)
      { exp_cov[l] = chk_alloc (N_train*N_train, sizeof(double));
        optional_left -= N_train*N_train*sizeof(double);
      }
      else
      { exp_cov[l] = 0;
      }
    }

    post_mean = chk_alloc (N_train, sizeof(double));

    scr1 = chk_alloc (N_train, sizeof(double));
    scr2 = chk_alloc (N_train, sizeof(double));

    /* Make sure we don't do all this again. */

    initialize_done = 1;
  }

  /* Set up Monte Carlo state structure. */

  ds->aux_dim = 2*gp->N_outputs*N_train;
  ds->aux     = aux;

  ds->dim = hypers.total_hypers;
  ds->q   = hypers.hyper_block;

  ds->temp_state = 0;
  
  ds->stepsize = stepsizes.hyper_block;
}


/* RESET INITIALIZE_DONE IN PREPARATION FOR NEW LOG FILE. */

void gp_mc_cleanup(void)
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
  logf->header.size = hypers.total_hypers * sizeof (double);
  log_file_append (logf, hypers.hyper_block);

  if (have_values)
  { logf->header.type = 'F';
    logf->header.index = index;
    logf->header.size = gp->N_outputs*N_train * sizeof (double);
    log_file_append (logf, latent_values);
  }

  if (have_variances)
  { logf->header.type = 'N';
    logf->header.index = index;
    logf->header.size = gp->N_outputs*N_train * sizeof (double);
    log_file_append (logf, noise_variances);
  }
}


/* APPLICATION-SPECIFIC SAMPLING PROCEDURE. */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double pm,
  mc_iter *it
)
{
  int i, j;

  if (strcmp(op,"scan-values")==0)
  {
    if (model==0)
    { fprintf (stderr,
       "The scan-values operation is not allowed when there is no model\n");
      exit(1);
    }

    if (!have_values)
    { for (i = 0; i<N_train; i++)
      { for (j = 0; j<gp->N_outputs; j++)
        { latent_values [gp->N_outputs*i + j] = 0;
        }
      }
      have_values = 1;
    }

    scan_values (pm==0 ? 1 : pm);

    ds->know_pot = 0;
    ds->know_grad = 0;

    return 1;
  }

  else if (strcmp(op,"sample-values")==0)
  {
    if (model==0)
    { fprintf (stderr,
       "The sample-values operation is not allowed when there is no model\n");
      exit(1);
    }

    if (model->type!='R')
    { fprintf (stderr,
       "The sample-values operation is not allowed for this model\n");
      exit(1);
    }

    sample_values(); 
    have_values = 1;

    ds->know_pot = 0;
    ds->know_grad = 0;

    return 1;
  }

  else if (strcmp(op,"discard-values")==0)
  {
    have_values = 0;

    ds->know_pot = 0;
    ds->know_grad = 0;

    return 1;
  }

  else if (strcmp(op,"sample-variances")==0)
  {
    if (model==0 || model->type!='R' || model->noise.alpha[2]==0)
    { fprintf (stderr,
       "The sample-variances operation is not allowed for this model\n");
      exit(1);
    }

    if (!have_values) 
    { sample_values(); /* Creates values, but doesn't remember them, since  */
    }                  /*   have_values remains set to zero.                */

    sample_variances();

    ds->know_pot = 0;
    ds->know_grad = 0;

    return 1;
  }

  return 0;
}


/* EVALUATE POTENTIAL ENERGY AND ITS GRADIENT.  Assumes that the dynamical
   state is the same as that in the 'hypers' variable of this module.  The
   energy is based on the prior for the hyperparameters and the likelihood.  
   The likelihood is based on the target values for training cases if no
   latent values are stored, and on the latent values for training cases if 
   these do exist.  Note that the latent values are considered fixed at this 
   point (they are changed by sample-values).  If case-by-case noise variances 
   are present, they too are considered fixed when computing the likelihood. */

void mc_app_energy
( mc_dynamic_state *ds,	/* Current dynamical state */
  int N_approx,		/* Number of gradient approximations in use */
  int w_approx,		/* Which approximation to use this time */
  double *energy,	/* Place to store energy, null if not required */
  mc_value *gr		/* Place to store gradient, null if not required */
)
{
  double ld, lp, found, c, a, n, v;
  int i, j, k;

  /* Set things up to start. */

  if (ds->q!=hypers.hyper_block) abort();

  if (gr && gr!=grad.hyper_block)
  { grad.total_hypers = hypers.total_hypers;
    grad.hyper_block = gr;
    gp_hyper_pointers (&grad, gp, model);
  }

  /* Evaluate portion of energy and gradient due to prior.  If the prior
     probability is very small, go to the recovery code. */

  lp = gp_log_prior(&hypers,gp,model,0);

  if (lp<-1e30) goto recover;

  if (energy) 
  { *energy = -lp;
  }

  if (gr)
  { gp_prior_grad(&hypers,gp,model,&grad);
  }

  /* Evaluate the portion of the energy and gradient due to the likelihood
     for the case-by-case noise variances, if present. */

  if (N_train>0 && have_variances 
       && model!=0 && model->type=='R' && model->noise.alpha[2]!=0)
  {
    a = model->noise.alpha[2];

    for (j = 0; j<gp->N_outputs; j++)
    { 
      n = exp(*hypers.noise[j]);

      for (i = 0; i<N_train; i++)
      { v = log(noise_variances[gp->N_outputs*i+j]) / 2;
        if (energy) 
        { *energy -= gp_gdens (a, n, v, 0);
        }
        if (gr && (model->noise.alpha[0]!=0 || model->noise.alpha[1]!=0)) 
        { gp_gdiff (a, n, v, grad.noise[j], 0);
        }
      }
    }
  }

  /* Evaluate portion of the energy and gradient due to the likelihood
     for target values give latent values, if present. */

  if (N_train>0 && have_values && !have_variances 
       && model!=0 && model->type=='R' && model->noise.alpha[2]==0)
  {
    for (j = 0; j<gp->N_outputs; j++)
    { 
      n = exp(*hypers.noise[j]);

      if (energy)
      { *energy += N_train*log(n);
      }

      if (gr && (model->noise.alpha[0]!=0 || model->noise.alpha[1]!=0))
      { *grad.noise[j] += N_train;
      }

      for (i = 0; i<N_train; i++)
      { v = latent_values[gp->N_outputs*i+j] - train_targets[gp->N_outputs*i+j];
        if (energy) 
        { *energy += (v*v)/(2*n*n);
        }
        if (gr && (model->noise.alpha[0]!=0 || model->noise.alpha[1]!=0)) 
        { *grad.noise[j] += -(v*v)/(n*n);
        }
      }
    }
  }

  /* Evaluate portion of energy and gradient due to the data likelihood.  
     This is done separately for each output, but with the computation of the
     covariance matrix done only once if there is no change from output
     to output. */

  if (N_train>0)
  {
    for (j = 0; j<gp->N_outputs; j++)
    { 
      /* Compute covariance matrix for this output, if it's not the same as 
         for the previous output, and then find its Cholesky decomposition
         (and log determinant), and its inverse.  If this process fails, 
         go to the recovery code. */

      if (j==0 || model!=0 && model->type=='R' && !have_values
        && (model->noise.alpha[2]!=0 || *hypers.noise[j]!=*hypers.noise[j-1]))
      {
        gp_cov (gp, &hypers, train_inputs, N_train, train_inputs, N_train,
                train_cov, gr ? exp_cov : 0);

        if (gp->has_jitter)
        { for (i = 0; i<N_train; i++)
          { train_cov[i*N_train+i] += exp(2 * *hypers.jitter);
          }
        }

        if (!have_values && model!=0)
        { if (model->type!='R') abort();
          if (model->noise.alpha[2]!=0)
          { if (!have_variances) abort();
            for (i = 0; i<N_train; i++) 
            { train_cov[i*N_train+i] += noise_variances[gp->N_outputs*i+j];
            }
          }
          else
          { for (i = 0; i<N_train; i++) 
            { train_cov[i*N_train+i] += exp(2 * *hypers.noise[j]);
            }
          }
        }

        if (!cholesky(train_cov,N_train,&ld)
         || !inverse_from_cholesky (train_cov, scr1, scr2, N_train))
        { goto recover;
        }
      }

      /* Multiply values for this output by inverse covariance; save in scr1. */

      for (i = 0; i<N_train; i++)
      { scr1[i] = inner_product(train_cov+i*N_train, 1, 
                     have_values ? latent_values+j : train_targets+j, 
                     gp->N_outputs, N_train);
      }

      /* Add the contribution of this output to the energy (if wanted).  Go
         to the recovery code if the contribution is huge. */

      c = ld/2 + inner_product(scr1, 1, 
                    have_values ? latent_values+j : train_targets+j, 
                    gp->N_outputs, N_train) / 2;

      if (c>1e30) goto recover;

      if (energy) *energy += c;

      /* Compute the contribution of this output to the gradient, if wanted.
         The derivative with respect to each hyperparameter is computed 
         separately. */

      if (gr)
      {
        for (k = 0; k<ds->dim; k++)
        { 
          /* Compute the derivative of the covariance matrix for output j with
             respect to hyperparameter k.  Initially, only the upper triangular
             part is computed; the lower part is filled in from this later. */

          found = gp_cov_deriv (gp, &hypers, exp_cov, ds->q+k, 
                                train_inputs, diff_cov, N_train);

          if (gp->has_jitter && ds->q+k==hypers.jitter)
          { for (i = 0; i<N_train; i++)
            { diff_cov[i*N_train+i] += 2 * exp(2 * *hypers.jitter);
            }
            found = 1;
          }

          if (model!=0 && model->type=='R' && !have_values)
          { if (ds->q+k==hypers.noise[j])
            { for (i = 0; i<N_train; i++) 
              { diff_cov[i*N_train+i] += 2 * exp(2 * *hypers.noise[j]);
              }
              found = 1;
            }
          }

          /* Add to derivative of energy using this matrix of derivatives,
             if this hyperparameter plays a direct role in the likelihood. */

          if (found)
          {
            fill_lower_triangle (diff_cov, N_train);

            matrix_product (diff_cov, scr1, scr2, N_train, 1, N_train);

            gr[k] += trace_of_product(train_cov, diff_cov, N_train) / 2
                   - inner_product(scr1, 1, scr2, 1, N_train) / 2;
          }
        }
      }
    }
  }

  /* Check for gradient being huge, and suppress if so. */

  if (gr)
  { for (k = 0; k<ds->dim; k++)
    { if (gr[k]>1e10 || gr[k]<-1e10) 
      { goto recover_gr;
      }
    }
  }

  return;

  /* Recovery code for when the state is apparently ridiculous.  Just 
     set the energy to a huge value (so we will reject), and set the gradient
     to zero (so we don't make things even worse).  The recover_gr entry
     is for when the situation was diagnosed based on the gradient (so for
     consistency, we shouldn't alter the computed energy). */

recover:

  if (energy) 
  { *energy = 1e30;
  }

recover_gr:

  if (gr)
  { for (k = 0; k<ds->dim; k++) gr[k] = 0;
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


/* SET STEPSIZES FOR EACH COORDINATE.  Assumes that the stepsizes are
   stored in the same place as the 'stepsizes' variable of this module. 
   Stepsizes are set heuristically based on the assumption that the
   range of hyperparameters other than the noise variances does not
   get smaller as the training set size increases. This is somewhat
   debatable. */

void mc_app_stepsizes
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
  int i, l;

  if (ds->stepsize!=stepsizes.hyper_block) abort();

  if (gp->has_constant)
  { 
    if (gp->constant.alpha[0]!=0)
    { *stepsizes.constant = 0.1;
    }
  }

  if (gp->has_linear)
  { 
    if (gp->linear.alpha[0]!=0)
    { *stepsizes.linear_cm = 0.1/sqrt(gp->N_inputs);
    }

    if (gp->linear.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { *stepsizes.linear[i] = 0.1;
      }
    }
  }

  if (gp->has_jitter)
  { 
    if (gp->jitter.alpha[0]!=0)
    { *stepsizes.jitter = 1/sqrt(N_train*gp->N_outputs+1);
    }
  }

  for (l = 0; l<gp->N_exp_parts; l++)
  { 
    if (gp->exp[l].scale.alpha[0]!=0)
    { *stepsizes.exp[l].scale = 0.1;
    }

    if (gp->exp[l].relevance.alpha[0]!=0)
    { *stepsizes.exp[l].rel_cm = 0.1/sqrt(gp->N_inputs);
    }

    if (gp->exp[l].relevance.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { *stepsizes.exp[l].rel[i] = 0.1;
      }
    }
  }

  if (model!=0 && model->type=='R')
  {  
    if (model->noise.alpha[0]!=0)
    { *stepsizes.noise_cm = 
         model->noise.alpha[1]==0 ? 1/sqrt(N_train*gp->N_outputs+1)
                                  : 0.1/sqrt(gp->N_outputs);
    }
    
    if (model->noise.alpha[1]!=0)
    { for (i = 0; i<gp->N_outputs; i++)
      { *stepsizes.noise[i] = 1/sqrt(N_train+1);
      }
    }
  }
}


/* LOG PROBABILITY FUNCTION USED IN IN SCAN_VALUES PROCEDURE.  Computes
   the log probability and the derivative of the log probability for
   a latent value that determines a binary or class probability, as 
   needed for the Adaptive Rejection Sampling procedure (ars).  Constant
   factors in the likelihood are ignored.  The distribution is described 
   by the extra structure below, which gives the prior mean and variance, 
   and the information needed to determine the likelihood. */

struct extra { int b; double cnst, mean, var; };

static double logp
( double x,		/* Point to evaluate log density at */
  double *d,		/* Place to store derivative of log density */
  void *extra		/* Extra information */
)
{
  struct extra *ex = extra;
  double e;

  e = exp(x);

  *d = - (x-ex->mean) / ex->var 
       - e / (ex->cnst + e)
       + (ex->b ? 1 : 0);

  return - (x-ex->mean)*(x-ex->mean) / (2*ex->var) 
         - log (ex->cnst + e) 
         + (ex->b ? x : 0);
}


/* DO GIBBS SAMPLING SCANS FOR LATENT VALUES. */

static void scan_values
( int scans	/* Number of Gibbs sampling scans to do for each output */
)
{
  double mean, prec;
  int i, j, k;
  int N_outputs;

  N_outputs = gp->N_outputs;

  /* Compute the covariance matrix for latent values at training points. 
     This will be the same for all outputs, since it doesn't include any
     noise part. */

  gp_cov (gp, &hypers, train_inputs, N_train, train_inputs, N_train,
          latent_cov, 0);

  if (gp->has_jitter)
  { for (i = 0; i<N_train; i++)
    { latent_cov[i*N_train+i] += exp(2 * *hypers.jitter);
    }
  }

  /* Invert this covariance matrix. */

  if (!cholesky(latent_cov,N_train,0)
   || !inverse_from_cholesky (latent_cov, scr1, scr2, N_train))
  { fprintf(stderr,"Couldn't find inverse of covariance in scan-values!\n");
    exit(1);
  }

  /* Do Gibbs sampling scans for latent values. */

  for (k = 0; k<scans; k++)
  {
    for (j = 0; j<N_outputs; j++)
    {
      for (i = 0; i<N_train; i++)
      {
        /* Find mean and precision (inverse variance) for this value, 
           based on other values. */

        latent_values[i*N_outputs+j] = 0;

        mean = - inner_product (latent_values+j, N_outputs, 
                   latent_cov+N_train*i, 1, N_train) / latent_cov[i*N_train+i];

        prec = latent_cov[i*N_train+i];

        /* Combine this with the likelihood from the target value and
           sample from the resulting distribution. */

        if (model->type=='R')
        { double pr;
          pr = have_variances ? 1/noise_variances[i*N_outputs+j]
                              : exp (- 2 * *hypers.noise[j]);
          mean = (mean*prec + train_targets[i*N_outputs+j]*pr) / (prec+pr);
          prec = prec+pr;
          latent_values[i*N_outputs+j] = mean + rand_gaussian()/sqrt(prec);
        }

        else if (model->type=='B')
        { struct extra ex;
          ex.cnst = 1;
          ex.mean = mean;
          ex.var  = 1/prec;
          ex.b    = train_targets[i*N_outputs+j]==1;
          latent_values[i*N_outputs+j] = ars (mean, sqrt(1/prec), logp, &ex);
        }

        else if (model->type=='C')
        { struct extra ex;
          int jj;
          ex.cnst = 0;
          for (jj = 0; jj<N_outputs; jj++)
          { if (jj!=j) ex.cnst += exp(latent_values[i*N_outputs+jj]);
          }
          ex.mean = mean;
          ex.var  = 1/prec;
          ex.b    = train_targets[i]==j;
          latent_values[i*N_outputs+j] = ars (mean, sqrt(1/prec), logp, &ex);
        }

        else
        { abort(); 
        }
      }
    }
  }
}


/* SAMPLE LATENT VALUES FOR A REGRESSION MODEL. */

static void sample_values(void)
{
  int N_outputs;
  int i, j;

  N_outputs = gp->N_outputs;

  for (j = 0; j<N_outputs; j++)
  { 
    /* Compute Cholesky decomposition of posterior covariance and 
       regression matrix for this output, if it's not the same as 
       for the last output. */

    if (j==0 || model->noise.alpha[2]!=0 
             || *hypers.noise[j]!=*hypers.noise[j-1])
    {
      /* Compute covariance matrix for this output, without added noise, and
         store in latent_cov. */

      gp_cov (gp, &hypers, train_inputs, N_train, train_inputs, N_train,
              latent_cov, 0);

      if (gp->has_jitter)
      { for (i = 0; i<N_train; i++)
        { latent_cov[i*N_train+i] += exp(2 * *hypers.jitter);
        }
      }

      /* Copy latent_cov to train_cov and add noise covariance. */

      for (i = N_train*N_train - 1; i>=0; i--) 
      { train_cov[i] = latent_cov[i];
      }

      if (model->noise.alpha[2]!=0)
      { if (!have_variances) abort();
        for (i = 0; i<N_train; i++) 
        { train_cov[i*N_train+i] += noise_variances[N_outputs*i+j];
        }
      }
      else
      { for (i = 0; i<N_train; i++) 
        { train_cov[i*N_train+i] += exp(2 * *hypers.noise[j]);
        }
      }
      
      /* Find Cholesky decomposition and use it to compute inverse. */
    
      if (!cholesky (train_cov,N_train,0)
       || !inverse_from_cholesky (train_cov, scr1, scr2, N_train))
      { fprintf (stderr, "Couldn't invert covariance in sample-values!\n");
        exit(1);
      }
    
      /* Find the matrix of coefficients for computing the posterior mean from
         the targets. */
    
      matrix_product (latent_cov, train_cov, reg_matrix, 
                      N_train, N_train, N_train);
    
      /* Find the posterior covariance for latent values underlying training 
         cases, using various matrices computed above. */
    
      matrix_product (reg_matrix, latent_cov, post_cov,
                      N_train, N_train, N_train);
    
      for (i = N_train*N_train - 1; i>=0; i--) 
      { post_cov[i] = latent_cov[i] - post_cov[i];
      }

      /* Find the Cholesky decomposition of the posterior covariance. */

      if (!cholesky(post_cov,N_train,0))
      { fprintf(stderr,
          "Couldn't find Cholesky decomposition of posterior covariance in sample-values!\n");
         exit(1);
      }
    }

    /* Sample latent values for each training case. */

    for (i = 0; i<N_train; i++) 
    { scr1[i] = rand_gaussian();
    }

    for (i = 0; i<N_train; i++)
    { latent_values[i*N_outputs+j] = 
         inner_product (scr1, 1, post_cov+i*N_train, 1, i+1)
          + inner_product (reg_matrix+i*N_train, 1, train_targets+j, 
                           N_outputs, N_train);
    }
  }
}


/* SAMPLE CASE-BY-CASE VARIANCES FOR REGRESSION MODEL. */

static void sample_variances(void)
{
  double a, n, d, s;
  int N_outputs;
  int i, j;

  N_outputs = gp->N_outputs;

  a = model->noise.alpha[2];

  for (j = 0; j<N_outputs; j++)
  { 
    n = exp (2 * *hypers.noise[j]);

    for (i = 0; i<N_train; i++)
    {
      d = latent_values[N_outputs*i+j] - train_targets[N_outputs*i+j];
      s = rand_gamma((a+1)/2) / ((a*n+d*d)/2);

      noise_variances[N_outputs*i+j] = 1/s;
    }
  }
}
