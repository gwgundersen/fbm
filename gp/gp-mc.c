/* GP-MC.C - Interface between Gaussian process and Markov chain modules. */

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


/* USE CHEAP ENERGY FUNCTION?  If set to 0, the constant terms in minus the
   log likelihood are included in the energy; if set to 1, they are not.
   This is relevant if the energy is going to be interpreted to give the
   likelihood. */

#define Cheap_energy 0		/* Normally set to 0 */


/* SAVE NON-LINEAR (EG, EXPONENTIAL) TERMS IN COVARIANCE?  Set to zero
   only in order to debug the code that handles the case when there's
   not enough memory to do this. */

#define Use_exp_cov 1		/* Normally set to 1 */


/* CONSTANTS INVOLVING PI. */

#ifndef M_PI
#define M_PI 3.14159265358979323846	/* Define pi, if not defined already */
#endif

#define Log2pi  1.83787706640934548356	/* Log(2*M_PI) */


/* ADJUSTABLE PARAMETER CONTROLLING OPTIONAL MEMORY USAGE.  This adjustable
   parameter sets the maximum amount of memory that will be used to store
   intermediate results (specifically, terms in the covariances), in order to 
   avoid recomputing them.  Its setting has no effect on the results, but 
   setting it higher may reduce computation time, at the cost of memory usage.*/

#define Max_optional_memory 20000000  /* Max. optional memory usage, in bytes */


/* GAUSSIAN PROCESS VARIABLES. */

static int initialize_done = 0;	/* Has this all been set up? */

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

static double *diag_inv;	/* Diagonal of inverse covariance (N_train) */

static double *scr1, *scr2;	/* Scratch vectors (N_train long) */


/* PROCEDURES. */

static void sample_values (void);
static void jitter_values (double, int);
static void sample_variances (void);
static void scan_values (int);
static void met_values (double, int, mc_iter *);
static void mh_values (double, int, mc_iter *);


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

    if (model!=0 && model->type=='V')
    { fprintf(stderr,"Can't handle survival models in gp-mc\n");
      exit(1);
    }
  
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
      if (!Use_exp_cov)
      { exp_cov[l] = 0;
      }
      else if (l==0)
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

    diag_inv = chk_alloc (N_train, sizeof(double));

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


/* APPLICATION-SPECIFIC SAMPLING PROCEDURES. */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double pm,
  double pm2,
  mc_iter *it,
  mc_temp_sched *sch
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

    if (model->type=='R' && model->autocorr)
    { fprintf(stderr,
       "The scan-values operation isn't implemented for models with autocorrelated noise\n");
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

    scan_values (pm<=0 ? 1 : pm);

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

  else if (strcmp(op,"jitter-values")==0)
  {
    if (model==0)
    { fprintf (stderr,
       "The jitter-values operation is not allowed when there is no model\n");
      exit(1);
    }

    if (!gp->has_jitter)
    { fprintf (stderr,
       "The jitter-values operation is not allowed when there is no jitter\n");
      exit(1);
    }

    if (pm<0 || pm>=1)
    { fprintf(stderr,"First parameter after jitter-values must be in [0,1)\n");
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

    jitter_values (pm, (int)pm2<=0 ? 1 : (int)pm2); 
    have_values = 1;

    ds->know_pot = 0;
    ds->know_grad = 0;

    return 1;
  }

  else if (strcmp(op,"met-values")==0)
  { 
    if (model==0)
    { fprintf (stderr,
       "The met-values operation is not allowed when there is no model\n");
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

    if (pm<0)
    { fprintf(stderr,"First parameter after met-values must be positive\n");
      exit(1);
    }

    met_values (pm==0 ? 1 : pm, (int)pm2<=0 ? 1 : (int)pm2, it);

    ds->know_pot = 0;
    ds->know_grad = 0;

    return 1;
  }

  else if (strcmp(op,"mh-values")==0)
  { 
    if (model==0)
    { fprintf (stderr,
       "The mh-values operation is not allowed when there is no model\n");
      exit(1);
    }

    if (pm<0 || pm>1)
    { fprintf(stderr,"First parameter after mh-values must be in (0,1]\n");
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

    mh_values (pm==0 ? 1 : pm, (int)pm2<=0 ? 1 : (int)pm2, it);

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
  double ld, lp, c, a, n, v, s;
  int found;
  double *d;
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
     for target values given latent values, if present.  This is relevant
     only when it affect the noise variance hyperparameter for a regression
     model. */

  if (N_train>0 && have_values && !have_variances 
       && model!=0 && model->type=='R' && model->noise.alpha[2]==0)
  {
    if (model->autocorr)
    { fprintf(stderr,"Latent values must be discarded before hyperparameter\n");
      fprintf(stderr,"   updates are done when the noise is autocorrelated.\n");
      exit(1);
    }

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
    if (energy && !Cheap_energy) 
    { *energy += 0.5 * gp->N_outputs * N_train * Log2pi;
    }

    for (j = 0; j<gp->N_outputs; j++)
    { 
      /* Compute covariance matrix for this output, if it's not the same as 
         for the previous output, and then find its Cholesky decomposition
         (and log determinant), and if needed, its inverse (stored in lower
         part of diff_cov and in diag_inv).  If this process fails, go to the 
         recovery code. */

      if (j==0 || model!=0 && model->type=='R' && !have_values
        && (model->noise.alpha[2]!=0 || *hypers.noise[j]!=*hypers.noise[j-1]))
      {
        gp_train_cov (gp, model, &hypers, j, noise_variances, 
                      have_values ? train_cov : 0,
                      have_values ? 0 : train_cov,
                      gr ? exp_cov : 0);

        if (!cholesky(train_cov,N_train,&ld)) 
        { goto recover;
        }
 
        if (gr)
        { for (i = N_train*N_train-1; i>=0; i--)
          { diff_cov[i] = train_cov[i];
          }
          if (!inverse_from_cholesky (diff_cov, scr1, scr2, N_train))
          { goto recover;
          }
          for (i = 0; i<N_train; i++)
          { diag_inv[i] = diff_cov[i*N_train+i];
          }
        }
      }

      /* Multiply values for this output by inverse of Cholesky decomposition.
         Put in scr1. */

      forward_solve (train_cov, scr1, 1, 
         have_values ? latent_values+j : train_targets+j, gp->N_outputs, 
         N_train);

      /* Add the contribution of this output to the energy (if wanted).  Go
         to the recovery code if the contribution is huge. */

      c = ld/2 + squared_norm (scr1, 1, N_train) / 2;

      if (c>1e30) goto recover;

      if (energy) *energy += c;

      /* Compute the contribution of this output to the gradient, if wanted.
         The derivative with respect to each hyperparameter is computed 
         separately. */

      if (gr)
      {
        /* Compute inverse of covariance times values and put in scr2. */

        backward_solve (train_cov, scr2, 1, scr1, 1, N_train);

        for (k = 0; k<ds->dim; k++)
        { 
          /* Compute the derivative of the covariance matrix for output j with
             respect to hyperparameter k.  Only the upper triangular part is 
             computed, and stored in diff_cov, allowing the part below the
             diagonal to still contain elements of the inverse covariance. */

          found = gp_cov_deriv (gp, &hypers, exp_cov, ds->q+k, 
                                train_inputs, diff_cov, N_train);

          if (gp->has_jitter && ds->q+k==hypers.jitter)
          { for (i = 0; i<N_train; i++)
            { diff_cov[i*N_train+i] += 2 * exp(2 * *hypers.jitter);
            }
            found = 2;
          }

          if (model!=0 && model->type=='R' && !have_values && !have_variances)
          { if (ds->q+k==hypers.noise[j])
            { double n;
              int lag;
              n = 2 * exp(2 * *hypers.noise[j]);
              for (i = 0; i<N_train; i++) 
              { diff_cov[i*N_train+i] += n;
              }
              if (model->autocorr)
              { for (i = 0; i<N_train; i++)
                { for (lag = 1; lag<=model->n_autocorr && i+lag<N_train; lag++) 
                  { diff_cov[i*N_train+(i+lag)] += n * model->acf[lag-1];
                  }
                }
                found = 1;
              }
              else
              { found = 2;
              }
            }
          }

          /* Add to derivative of energy using this matrix of derivatives,
             if this hyperparameter plays a direct role in the likelihood. */

          if (found)
          {
            d = diff_cov;
            v = d[0] * (diag_inv[0] - scr2[0]*scr2[0]); 

            if (found==2) /* Non-zero derivatives on diagonal only */
            { for (i = 1; i<N_train; i++)
              { d += N_train;
                s = scr2[i];
                v += d[i] * (diag_inv[i] - s*s);
              }
            }
            else /* General case */
            { for (i = 1; i<N_train; i++)
              { d += N_train;
                s = scr2[i];
                v += d[i] * (diag_inv[i] - s*s)
                      + 2 * inner_product (d, 1, diff_cov+i, N_train, i)
                      - 2 * s * inner_product (scr2, 1, diff_cov+i, N_train, i);
              }
            }

            gr[k] += v / 2;
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


/* SAMPLE FROM DISTRIBUTION AT INVERSE TEMPERATURE OF ZERO.  Returns zero
   if this is not possible. */

int mc_app_zero_gen
( mc_dynamic_state *ds	/* Current dynamical state */
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
    { if (gp->linear.alpha[1]==0)
      { *stepsizes.linear_cm = 0.1;
      }
      else
      { int c;
        c = 0;
        for (i = 0; i<gp->N_inputs; i++) 
        { c += !(gp->linear_flags[i]&Flag_omit);
        }
        *stepsizes.linear_cm = 0.1/sqrt((double)c);
      }
    }
   
    if (gp->linear.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->linear_flags[i]&Flag_omit))
        { *stepsizes.linear[i] = 0.1;
        }
      }
    }
  }

  if (gp->has_jitter)
  { 
    if (gp->jitter.alpha[0]!=0)
    { *stepsizes.jitter = 0.5/sqrt((double)N_train*gp->N_outputs+1);
    }
  }

  for (l = 0; l<gp->N_exp_parts; l++)
  { 
    if (gp->exp[l].scale.alpha[0]!=0)
    { *stepsizes.exp[l].scale = 0.1;
    }

    if (gp->exp[l].relevance.alpha[0]!=0)
    { if (gp->exp[l].relevance.alpha[1]==0)
      { *stepsizes.exp[l].rel_cm = 0.1;
      }
      else
      { int c;
        c = 0;
        for (i = 0; i<gp->N_inputs; i++) 
        { c += !(gp->exp[l].flags[i]&Flag_omit);
        }
        *stepsizes.exp[l].rel_cm = 0.1/sqrt((double)c);
      }
    }

    if (gp->exp[l].relevance.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->exp[l].flags[i]&Flag_omit))
        { *stepsizes.exp[l].rel[i] = 0.1;
        }
      }
    }
  }

  if (model!=0 && model->type=='R')
  {  
    if (model->noise.alpha[0]!=0)
    { *stepsizes.noise_cm = 
         model->noise.alpha[1]==0 ? 0.5/sqrt((double)N_train*gp->N_outputs+1)
                                  : 0.1/sqrt((double)gp->N_outputs);
    }
    
    if (model->noise.alpha[1]!=0)
    { for (i = 0; i<gp->N_outputs; i++)
      { *stepsizes.noise[i] = 0.5/sqrt((double)N_train+1);
      }
    }
  }
}


/* LOG PROBABILITY FUNCTION USED IN SCAN-VALUES AND JITTER-VALUES.  Computes
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


/* ANOTHER LOG PROBABILITY FUNCTION USED IN SCAN-VALUES AND JITTER-VALUES.  
   Like logp above, but for Poisson models.  Negative values are left
   censored up to their absolute value. */

static double logpp
( double x,		/* Point to evaluate log density at */
  double *d,		/* Place to store derivative of log density */
  void *extra		/* Extra information */
)
{
  struct extra *ex = extra;
  double e, v, s, lg;
  int i;

  e = exp(x);

  *d = - (x-ex->mean) / ex->var;
  v = - (x-ex->mean)*(x-ex->mean) / (2*ex->var);

  if (ex->b>0)
  { *d += ex->b - e;
    v += ex->b*x - e;
  }
  else
  {
    s = -e;
    for (i = 1; i <= -ex->b; i++)
    { s = addlogs (s, i*x - e - lgamma(i+1));
    }
   
    v += s;

    s = 0;
    for (i = 0; i <= -ex->b; i++)
    { s += exp (lgamma(-ex->b+1) - lgamma(i+1) - (-ex->b+1-i)*x);
    }

    *d += -1/s;
  }

  return v;
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

  gp_train_cov (gp, 0, &hypers, 0, 0, latent_cov, 0, 0);

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

        else if (model->type=='N')
        { struct extra ex;
          ex.mean = mean;
          ex.var  = 1/prec;
          ex.b    = train_targets[i*N_outputs+j];
          latent_values[i*N_outputs+j] = ars (mean, sqrt(1/prec), logpp, &ex);
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


/* DO METROPOLIS UPDATES FOR LATENT VALUES.  The Metropolis proposal is
   a change by an amount that is a scaled down sample from the prior,
   except when 'scale' is negative, in which case the proposed changes to
   different values are independent, with each having a standard deviation 
   that is scaled down from its prior standard deviation.  When there is 
   more than one output, updates are proposed separately for each in turn. */

static void met_values
( double scale,		/* Amount by which to scale covariance for proposal */
  int repeat,		/* Number of repetitions to do */
  mc_iter *it
)
{  
  int i, j, k;
  int N_outputs;

  double new_E, old_E;

  N_outputs = gp->N_outputs;

  /* Compute the covariance matrix for latent values at training points. 
     This will be the same for all outputs. */

  gp_train_cov (gp, 0, &hypers, 0, 0, latent_cov, 0, 0);

  /* Find Cholesky decomposition for use in evaluating energy, and in 
     generating correlated proposals. */

  if (!cholesky(latent_cov,N_train,0))
  { fprintf(stderr,"Couldn't find Cholesky decomposition in met-values!\n");
    exit(1);
  }

  /* Do Metropolis updates. */

  for (k = 0; k<repeat; k++)
  {
    for (j = 0; j<N_outputs; j++)
    { 
      /* Find initial "energy", from prior and likelihood. */

      forward_solve (latent_cov, scr1, 1, latent_values+j, gp->N_outputs, 
                     N_train);
  
      old_E = squared_norm (scr1, 1, N_train) / 2;

      for (i = 0; i<N_train; i++)
      { old_E -= gp_likelihood (&hypers, model, data_spec, 
                  train_targets+i*data_spec->N_targets, 
                  latent_values+i*gp->N_outputs, 
                  have_variances ? noise_variances+i*data_spec->N_targets : 0);
      }

      /* Save current latent variables in scr2. */

      for (i = 0; i<N_train; i++) 
      { scr2[i] = latent_values[i*N_outputs+j];
      }

      /* Add random change to latent variables. */

      for (i = 0; i<N_train; i++) 
      { scr1[i] = rand_gaussian();
      }

      if (scale>0)
      { for (i = 0; i<N_train; i++)
        { latent_values[i*N_outputs+j] += 
            scale * inner_product (scr1, 1, latent_cov+i*N_train, 1, i+1);
        }
      }
      else
      { for (i = 0; i<N_train; i++)
        { latent_values[i*N_outputs+j] += 
            scale * sqrt(latent_cov[i*N_train+i]) * scr1[i];
        }
      }

      /* Find new "energy", from prior and likelihood. */

      forward_solve (latent_cov, scr1, 1, latent_values+j, gp->N_outputs, 
                     N_train);

      new_E = squared_norm (scr1, 1, N_train) / 2;

      for (i = 0; i<N_train; i++)
      { new_E -= gp_likelihood (&hypers, model, data_spec, 
                  train_targets+i*data_spec->N_targets, 
                  latent_values+i*gp->N_outputs, 
                  have_variances ? noise_variances+i*data_spec->N_targets : 0);
      }   

      /* Accept or reject proposed change based on delta. */

      it->proposals += 1;
      it->delta = new_E - old_E;

      if (it->delta<=0 || rand_uniform() < exp(-it->delta))
      { it->move_point = 1;
      }
      else
      { for (i = 0; i<N_train; i++)
        { latent_values[i*N_outputs+j] = scr2[i];
        }
        it->rejects += 1; 
        it->move_point = 0;
      }
    }
  }
}


/* DO METROPOLIS-HASTINGS UPDATES FOR LATENT VALUES.  The proposal is to 
   latent values found by multiplying the current values by a constant factor
   and then adding a random amount that is a scaled down sample from the prior.
   The multiplicative factor is chosen so that the prior is invariant under 
   this update.  This proposal is then accepted or rejected based on the 
   likelihoods alone.  When there is more than one output, updates are 
   proposed separately for each in turn. */

static void mh_values
( double scale,		/* Amount by which to scale random change */
  int repeat,		/* Number of repetitions to do */
  mc_iter *it
)
{  
  int i, j, k;
  int N_outputs;

  double new_E, old_E;
  double alpha;

  N_outputs = gp->N_outputs;
 
  alpha = sqrt(1-scale*scale);

  /* Compute the covariance matrix for latent values at training points. 
     This will be the same for all outputs. */

  gp_train_cov (gp, 0, &hypers, 0, 0, latent_cov, 0, 0);

  /* Find Cholesky decomposition for use in generating correlated noise. */

  if (!cholesky(latent_cov,N_train,0))
  { fprintf(stderr,"Couldn't find Cholesky decomposition in mh-values!\n");
    exit(1);
  }

  /* Do Metropolis-Hastings updates. */

  for (k = 0; k<repeat; k++)
  {
    for (j = 0; j<N_outputs; j++)
    { 
      /* Find initial "energy" from likelihood. */

      old_E = 0;

      for (i = 0; i<N_train; i++)
      { old_E -= gp_likelihood (&hypers, model, data_spec, 
                  train_targets+i*data_spec->N_targets, 
                  latent_values+i*gp->N_outputs, 
                  have_variances ? noise_variances+i*data_spec->N_targets : 0);
      }

      /* Save current latent variables in scr2. */

      for (i = 0; i<N_train; i++) 
      { scr2[i] = latent_values[i*N_outputs+j];
      }

      /* Multiply and add random change to latent variables. */

      for (i = 0; i<N_train; i++) 
      { scr1[i] = rand_gaussian();
      }

      for (i = 0; i<N_train; i++)
      { latent_values[i*N_outputs+j] = alpha * latent_values[i*N_outputs+j]
              + scale * inner_product (scr1, 1, latent_cov+i*N_train, 1, i+1);
      }

      /* Find new "energy" from likelihood. */

      new_E = 0;

      for (i = 0; i<N_train; i++)
      { new_E -= gp_likelihood (&hypers, model, data_spec, 
                  train_targets+i*data_spec->N_targets, 
                  latent_values+i*gp->N_outputs, 
                  have_variances ? noise_variances+i*data_spec->N_targets : 0);
      }   

      /* Accept or reject proposed change based on delta. */

      it->proposals += 1;
      it->delta = new_E - old_E;

      if (it->delta<=0 || rand_uniform() < exp(-it->delta))
      { it->move_point = 1;
      }
      else
      { for (i = 0; i<N_train; i++)
        { latent_values[i*N_outputs+j] = scr2[i];
        }
        it->rejects += 1; 
        it->move_point = 0;
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
      /* Compute covariance matrix for both latent and target values. */

      gp_train_cov (gp, model, &hypers, j, noise_variances, 
                    latent_cov, train_cov, 0);
      
      /* Find Cholesky decomposition of train_cov and compute its inverse. */
    
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


/* UPDATE LATENT VALUES BY PLAYING WITH JITTER-REDUCED VERSIONS. */

static void jitter_values
( double jp,
  int repeat
)
{
  double jitter, mean, prec;
  int N_outputs;
  int i, j, k;

  N_outputs = gp->N_outputs;

  jitter = exp(2 * *hypers.jitter);

  /* Compute the covariance matrix for full latent values at training points, 
     and store in latent_cov.  This will be the same for all outputs, since 
     it doesn't include any noise part. */

  gp_train_cov (gp, 0, &hypers, 0, 0, latent_cov, 0, 0);

  /* Find covariance matrix of jitter-reduced latent variables, and store
     in train_cov. */

  for (i = N_train*N_train - 1; i>=0; i--) 
  { train_cov[i] = latent_cov[i];
  }

  for (i = 0; i<N_train; i++)
  { train_cov[i*N_train+i] -= (1-jp) * jitter;
  }

  /* Find Cholesky decomposition of the covariance of the full latent 
     variables and compute its inverse. */
    
  if (!cholesky (latent_cov,N_train,0)
   || !inverse_from_cholesky (latent_cov, scr1, scr2, N_train))
  { fprintf (stderr, "Couldn't invert covariance in jitter-values!\n");
    exit(1);
  }

  /* Find the matrix of coefficients for computing the mean of the 
     jitter-reduced latent variables from the full latent variables. */
    
  matrix_product (train_cov, latent_cov, reg_matrix, N_train, N_train, N_train);
    
  /* Find the covariance for jitter-reduced latent values conditional on
     the full latent variables, using various matrices computed above. */
    
  matrix_product (reg_matrix, train_cov, post_cov, N_train, N_train, N_train);
    
  for (i = N_train*N_train - 1; i>=0; i--) 
  { post_cov[i] = train_cov[i] - post_cov[i];
  }

  /* Find the Cholesky decomposition of this covariance, to use in generation.*/

  if (!cholesky(post_cov,N_train,0))
  { fprintf(stderr,
      "Couldn't find Cholesky decomposition of covariance in jitter-values!\n");
     exit(1);
  }

  for (k = 0; k<repeat; k++)
  {
    for (j = 0; j<N_outputs; j++)
    { 
      /* Sample jitter-reduced latent values for each training case. */

      for (i = 0; i<N_train; i++) 
      { scr1[i] = rand_gaussian();
      }

      for (i = 0; i<N_train; i++)
      { scr2[i] = inner_product (scr1, 1, post_cov+i*N_train, 1, i+1)
                + inner_product (reg_matrix+i*N_train, 1, latent_values+j, 
                                 N_outputs, N_train);
      }

      /* Sample full latent values based on jitter-reduced latent values
         and targets. */

      for (i = 0; i<N_train; i++)
      {
        /* Find mean and precision (inverse variance) for full latent value
           based on jitter-reduced value. */

        mean = scr2[i];
        prec = 1 / ((1-jp)*jitter);

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

        else if (model->type=='N')
        { struct extra ex;
          ex.mean = mean;
          ex.var  = 1/prec;
          ex.b    = train_targets[i*N_outputs+j];
          latent_values[i*N_outputs+j] = ars (mean, sqrt(1/prec), logpp, &ex);
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
