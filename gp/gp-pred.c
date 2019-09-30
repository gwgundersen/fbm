/* GP-PRED.C - Make predictions for test cases using Gaussian processes. */

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
 *
 * Features allowing selection of a given number of iterations and 
 * specification of ranges by cpu-time are adapted from code written
 * by Carl Edward Rasmussen, 1995.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "misc.h"
#include "rand.h"
#include "matrix.h"
#include "log.h"
#include "prior.h"
#include "model.h"
#include "data.h"
#include "numin.h"
#include "gp.h"
#include "gp-data.h"
#include "mc.h"
#include "pred.h"


/* NAME OF THIS MODULE. */

char *pred_app_name = "gp";


/* LOCAL VARIABLES. */

#define Prediction_sample 100	/* Size of sample for finding probabilities
				   for binary and class models */
static gp_spec *gp;

static gp_hypers hypers, *h = &hypers;
 
static int have_values, have_variances;

static double *latent_values, *noise_variances;

static double *train_cov, *scr1, *scr2, *scr3;

static double guessp, targ, prb, *prb_dist1, *prb_dist2;
static double *meanp, *varp; 


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* GAUSSIAN CUMULATIVE DISTRIBUTION FUNCTION. */

#define Phi(z) (0.5*erf((z)/sqrt(2.0))+0.5)


/* SET SIZES FOR APPLICATION RECORDS. */

void pred_app_record_sizes (void)
{
  gp_record_sizes(&logg);
}


/* INITIALIZE APPLICATION PROCEDURES. */

void pred_app_init (void)
{
  if (op_p && m==0)
  { fprintf(stderr,"Illegal combination of options with data model\n");
    exit(1);
  } 
 
  if (op_N && keep[0])
  { fprintf(stderr,
"Option '0' not allowed: Non-linear components are numbered starting with '1'\n"
     );
    exit(1);
  }

  if (m!=0 && m->type=='V')
  { fprintf(stderr,"Can't handle survival models in gp-pred\n");
    exit(1);
  }

  gp = logg.data['P'];

  gp_check_specs_present(gp,0,m,sv);

  h->total_hypers = gp_hyper_count(gp,m);

  logg.req_size['S'] = h->total_hypers * sizeof(double);

  gp_data_read (1, 1, gp, m, sv);

  train_cov = chk_alloc (N_train*N_train, sizeof(double));

  scr1 = chk_alloc (N_train, sizeof(double));
  scr2 = chk_alloc (N_train, sizeof(double));
  scr3 = chk_alloc (N_train, sizeof(double));

  meanp = chk_alloc (N_test*gp->N_outputs, sizeof(double));
  varp  = chk_alloc (N_test*gp->N_outputs, sizeof(double));

  prb_dist1    = chk_alloc (M_targets, sizeof (double));
  prb_dist2    = chk_alloc (M_targets, sizeof (double));
}


/* LOOK AT GP STORED AT THE CURRENT INDEX. Returns 1 if there really
   something here, zero if not. */

int pred_app_use_index (void)
{    
  int want_variance;
  int i, j, k;

  /* See if there's really something here. */

  if (logg.index['S']!=logg.last_index)
  { 
    return 0;
  }

  /* See what we have sitting in the log file for this iteration. */

  h->hyper_block = logg.data['S'];

  gp_hyper_pointers (h, gp, m);

  have_values    = logg.data['F']!=0 
                     && logg.index['F']==logg.last_index;
  have_variances = logg.data['N']!=0 
                     && logg.index['N']==logg.last_index;

  if (m!=0 && m->type=='R' && !op_l) 
  { have_values = 0;
  }

  latent_values   = logg.data['F'];
  noise_variances = logg.data['N'];

  if (have_values
   && logg.actual_size['F']!=N_train*gp->N_outputs*sizeof(double))
  { fprintf(stderr,"Record with latent values is wrong length\n");
    exit(1);
  }

  if (have_variances
   && logg.actual_size['N']!=N_train*gp->N_outputs*sizeof(double))
  { fprintf(stderr,"Record with noise variances is wrong length\n");
    exit(1);
  }

  /* Find predictive means and (maybe) variances for each output
     of every test case, storing them in meanp and varp. */

  for (j = 0; j<gp->N_outputs; j++)
  {
    /* See if we want the variance for this output. */

    want_variance = op_p || op_d || op_q 
                     || op_r && data_spec->trans[data_spec->N_inputs+j].take_log
                     || m!=0 && m->type!='R';

    /* Compute covariance matrix for training cases, and its Cholesky
       decomposition and maybe from that its inverse. Don't redo computations 
       if previous matrix is still applicable.*/

    if (j==0 || m!=0 && m->type=='R' && !have_values
      && (m->noise.alpha[2]!=0 || *hypers.noise[j]!=*hypers.noise[j-1]))
    {
      gp_train_cov (gp, m, h, j, noise_variances, 
                    have_values ? train_cov : 0,
                    have_values ? 0 : train_cov,
                    0);

      if (!cholesky (train_cov, N_train, 0))
      { fprintf(stderr,
"Couldn't find Cholesky decomposition of covariance matrix for training cases\n"
        );
        exit(1);
      }

      if (use_inverse)
      { if (!inverse_from_cholesky (train_cov, scr1, scr2, N_train))
        { fprintf (stderr,
           "Couldn't invert covariance matrix of training cases\n");
          exit(1);
        }
      }
    }

    /* Pre-compute product of inverse covariance for training points
       with training targets or latent values, for use in finding 
       predictive means.  If we're doing the alternate mean computation,
       we pre-compute the product of the inverse of the Cholesky decomposition
       with the training targets.  Put in scr1. */

    if (alt_mean)
    { forward_solve (train_cov, scr1, 1, 
              have_values ? latent_values+j : train_targets+j, gp->N_outputs, 
              N_train);
    }
    else if (use_inverse)
    { for (i = 0; i<N_train; i++)
      { scr1[i] = inner_product (train_cov+i*N_train, 1, 
                     have_values ? latent_values+j : train_targets+j, 
                     gp->N_outputs, N_train);
      }
    }
    else 
    { forward_solve (train_cov, scr2, 1, 
              have_values ? latent_values+j : train_targets+j, gp->N_outputs, 
              N_train);
      backward_solve (train_cov, scr1, 1, scr2, 1, N_train);
    }

    /* Find the predictive mean/variance for this output for each 
       test case. */

    for (i = 0; i<N_test; i++)
    {
      /* Find covariances between training cases and test case.  Put in scr2. */

      gp_cov (gp, h, test_inputs+gp->N_inputs*i, 1, 
              train_inputs, N_train, scr2, 0, op_N ? keep+1 : 0);

      /* If needed, find product of these covariances with inverse of Cholesky
         decomposition, unless we're using the actual inverse.  Put in scr3. */

      if (!use_inverse && (want_variance || alt_mean))
      { forward_solve (train_cov, scr3, 1, scr2, 1, N_train);
      }

      /* Find predictive mean using the vectors precomputed above. */

      if (alt_mean)
      { meanp[i*gp->N_outputs+j] = inner_product (scr1, 1, scr3, 1, N_train);
      }
      else
      { meanp[i*gp->N_outputs+j] = inner_product (scr1, 1, scr2, 1, N_train);
      }

      /* Find the predictive variance if it will be needed. */

      if (want_variance)
      { 
        gp_cov (gp, h, test_inputs+gp->N_inputs*i, 1,
                       test_inputs+gp->N_inputs*i, 1, 
                       &varp[i*gp->N_outputs+j], 0, op_N ? keep+1 : 0);

        if (gp->has_jitter)
        { varp[i*gp->N_outputs+j] += exp(2 * *h->jitter);
        }

        if (use_inverse)
        { matrix_product (scr2, train_cov, scr3, 1, N_train, N_train);
          varp[i*gp->N_outputs+j] -= inner_product (scr3, 1, scr2, 1, N_train);
        }
        else
        { varp[i*gp->N_outputs+j] -= squared_norm (scr3, 1, N_train);
        }

        if (varp[i*gp->N_outputs+j]<=0)
        { fprintf(stderr,
   "WARNING: Predicted variance not positive (case %d, output %d, var %.1le)\n",
           i, j, varp[i*gp->N_outputs+j]);
          varp[i*gp->N_outputs+j] = 1e-30;
          fprintf(stderr,"          - replaced by 1e-30\n");
        }
      }
    }
  }

  /* Use the predictive means and variances to make predictions
     for each test case. */

  for (i = 0; i<N_test; i++) 
  { 
    if (m==0 || m->type=='R') /* Model for real data */
    {
      test_log_prob[i] = 0;

      for (j = 0; j<gp->N_outputs; j++)
      { 
        if (m!=0)
        { if (m->noise.alpha[2]!=0)
          { double n;
            n = prior_pick_sigma(exp(*h->noise[j]),m->noise.alpha[2]);
            varp[i*gp->N_outputs+j] += n*n;
          }
          else
          { varp[i*gp->N_outputs+j] += exp(2 * *h->noise[j]);
          }
        }

        if (op_r) 
        { 
          tr = &data_spec->trans[data_spec->N_inputs+j];

          meanp[i*gp->N_outputs+j] = 
            data_inv_trans(meanp[i*gp->N_outputs+j], *tr);

          if (tr->take_log)
          {
            if (op_n && m!=0 && m->type=='R' && (m->noise.alpha[0]!=0 
                 || m->noise.alpha[1]!=0 || m->noise.alpha[2]!=0))
            { fprintf(stderr,"Predictive mean is undefined\n");
              exit(1);
            }
  
            meanp[i*gp->N_outputs+j] 
              *= exp (varp[i*gp->N_outputs+j]/(tr->scale*tr->scale*2));
          }
        }

        test_targ_pred[M_targets*i+j] = meanp[i*gp->N_outputs+j];

        if (op_D)
        { test_targ_med[M_targets*i+j] = meanp[i*gp->N_outputs+j];
        }
      
        if (have_targets && op_p) 
        { 
          targ = test_targets[gp->N_outputs*i+j];

          test_log_prob[i] += - log(2*M_PI*varp[i*gp->N_outputs+j])/2 
            - ((targ-meanp[i*gp->N_outputs+j])
               * (targ-meanp[i*gp->N_outputs+j])) 
                 / (2*varp[i*gp->N_outputs+j]);

          if (op_r)
          { tr = &data_spec->trans[data_spec->N_inputs+j];
            test_log_prob[i] += log(tr->scale);
            if (tr->take_log)
            { test_log_prob[i] -= log (data_inv_trans 
                (test_targets[data_spec->N_targets*i+j], *tr));
            }
          }
        }

        if (op_d || op_q)
        {
          for (k = 0; k<Median_sample; k++)
          {
            guessp = meanp[i*gp->N_outputs+j] 
                      + sqrt(varp[i*gp->N_outputs+j])*rand_gaussian();

            if (op_r)
            { guessp = data_inv_trans (guessp, 
                                       data_spec->trans[data_spec->N_inputs+j]);
            }

            median_sample[i][j][ms_count+k] = guessp;
          }
        }
      }
    }

    else if (m->type=='B') /* Model for binary data */
    {
      test_log_prob[i] = 0;
 
      for (j = 0; j<gp->N_outputs; j++)
      {
        double mean, sd, av0, av1, pr0;
        int n0, n1;

        mean = meanp[i*gp->N_outputs+j];
        sd   = sqrt(varp[i*gp->N_outputs+j]);

        av0 = 0;
        av1 = 0.1;

        n0 = n1 = 0;

        for (k = 0; k<Prediction_sample; k++)
        { double v;
          v = mean + rand_gaussian() * sd;
          if (v<0) { n0 += 1; av0 += 1/(1+exp(-v)); }
          else     { n1 += 1; av1 += 1/(1+exp(-v)); }
        }

        av0 /= n0+0.1;
        av1 /= n1+0.1;

        pr0 = Phi(-mean/sd);

        prb = pr0*av0 + (1-pr0)*av1;

        test_targ_pred[M_targets*i+j] = prb;

        if (have_targets && op_p)
        { targ = test_targets[gp->N_outputs*i+j];
          test_log_prob[i] += targ==1 ? log(prb) : log(1-prb);
        }
      }
    }

    else if (m->type=='N') /* Model for count data */
    {
      test_log_prob[i] = 0;
 
      for (j = 0; j<gp->N_outputs; j++)
      {
        double mean, sd, tfact, prb;
        int targ;

        mean = meanp[i*gp->N_outputs+j];
        sd   = sqrt(varp[i*gp->N_outputs+j]);

        if (have_targets && op_p)
        { targ = test_targets[gp->N_outputs*i+j];
          tfact = lgamma(targ+1);
          prb = 0;
        }
        else
        { targ = -1;
        }

        test_targ_pred[M_targets*i+j] = 0;

        if (op_R)
        { test_targ_pred[M_targets*i+j] = mean;
        }

        if (!op_R || targ>=0)
        {
          for (k = 0; k<Prediction_sample; k++)
          { double v;
            v = mean + rand_gaussian() * sd;
            if (!op_R) 
            { test_targ_pred[M_targets*i+j] += exp(v);
            }
            if (targ>=0)
            { prb += exp (targ*v - tfact - exp(v));
            }
          }
    
          if (!op_R) 
          { test_targ_pred[M_targets*i+j] /= Prediction_sample;
          }
          if (targ>=0) 
          { test_log_prob[i] += log(prb/Prediction_sample);
          }
        }
      }
    }

    else if (m->type=='C') /* Model for multi-class data */
    {
      test_log_prob[i] = 0;
 
      for (j = 0; j<gp->N_outputs; j++) 
      { prb_dist1[j] = 0;
      }

      for (k = 0; k<Prediction_sample; k++)
      { double v, s;
        s = 0;
        for (j = 0; j<gp->N_outputs; j++)
        { v = exp (meanp[i*gp->N_outputs+j] 
                    + rand_gaussian() * sqrt(varp[i*gp->N_outputs+j]));
          prb_dist2[j] = v;
          s += v;
        }
        for (j = 0; j<gp->N_outputs; j++)
        { prb_dist2[j] /= s;
          prb_dist1[j] += prb_dist2[j];
        }
      }

      for (j = 0; j<gp->N_outputs; j++) 
      { prb_dist1[j] /= Prediction_sample;
      }

      for (j = 0; j<M_targets; j++) 
      { test_targ_pred[M_targets*i+j] = prb_dist1[j];
      }

      if (have_targets && op_p)
      { targ = test_targets[i];
        test_log_prob[i] += log(prb_dist1[(int)targ]);
      }
    }

    else /* Type of model that isn't handled */
    { 
      abort();
    }
  }

  return 1;
}


/* CLEAN UP WHEN END OF LOG FILE IS REACHED. */

void pred_app_finish_file (void)
{
  logg.index['S'] = -1;  /* So they won't be mistaken for records from */
  logg.index['V'] = -1;  /* the new log file.                          */
  logg.index['F'] = -1;
}
