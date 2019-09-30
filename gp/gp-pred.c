/* GP-PRED.C - Program to predict for test cases using Gaussian processes. */

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


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* SOME CONSTANTS. */

#define Max_median_points 200	/* Max points that can be used for median */
#define Median_sample 11	/* Size of median sample per process */

#define Prediction_sample 100	/* Size of sample for finding probabilities
				   for binary and class models */

/* GAUSSIAN CUMULATIVE DISTRIBUTION FUNCTION. */

#define Phi(z) (0.5*erf((z)/sqrt(2.0))+0.5)


/* LOCAL PROCEDURES. */

static float find_median (float *, int);

static void usage(void);


/* MAIN PROGRAM. */

void main
( int argc,
  char **argv
)
{
  gp_spec *gp;
  model_specification *m;
  model_survival *sv;

  gp_hypers hypers, *h = &hypers;

  mc_iter *it;

  int op_i, op_t, op_r, op_p, op_m, op_n, op_d, op_l, op_b, op_a;
  int put_in_blanks;
  int guessing;
  char *options;

  int M_targets, N_gps;

  double *sum_targets;

  double *sq_error, *sq_error_sq, *abs_error, *abs_error_sq, *wrong, *wrong_sq;
  double tsq_error, tsq_error_sq, tabs_error, tabs_error_sq;
  double *guess, *error;

  double *log_prob, ave_log_prob, ave_log_prob_sq;

  float ***median_sample;
  int ms_count;

  log_file logf;
  log_gobbled logg;

  char test_inputs_spec[100];
  char test_targets_spec[100];
 
  int have_values, have_variances;
  int have_targets;

  double *latent_values, *noise_variances;

  double *train_cov, *scr1, *scr2, *scr3;

  double curr_log_prob, guessp, targ, prb, *prb_dist1, *prb_dist2;
  double *meanp, *varp; 

  data_transformation *tr;

  double dlindex, dhindex, target_index;
  int lindex, hindex, mod_no;
  char **ap;
  int i, j, k;

  /* Look at arguments. */

  if (argc<3) usage();

  test_inputs_spec[0] = 0;
  test_targets_spec[0] = 0;

  if (strcmp(argv[argc-2],"/")==0)
  { strcpy(test_inputs_spec,argv[argc-1]);
    argc -= 2;
  }
  else if (strcmp(argv[argc-3],"/")==0)
  { strcpy(test_inputs_spec,argv[argc-2]);
    strcpy(test_targets_spec,argv[argc-1]);
    argc -= 3;
  }

  argv[argc] = 0;

  if (argc<3) usage();

  options = argv[1];

  op_i = strchr(options,'i')!=0;
  op_t = strchr(options,'t')!=0;
  op_r = strchr(options,'r')!=0;
  op_p = strchr(options,'p')!=0;
  op_m = strchr(options,'m')!=0;
  op_n = strchr(options,'n')!=0;
  op_d = strchr(options,'d')!=0;
  op_l = strchr(options,'l')!=0;
  op_b = strchr(options,'b')!=0 || strchr(options,'B')!=0;
  op_a = strchr(options,'a')!=0;

  if (strlen(options) != op_i+op_t+op_r+op_p+op_m+op_n+op_d+op_l+op_b+op_a)
  { usage();
  }

  guessing = op_p || op_m || op_n || op_d;

  put_in_blanks = strchr(options,'B')!=0;

  if ((argc-2)%2!=0 && (guessing || argc!=3)) usage();

  ap = argv+2;

  /* Open first log file and read Gaussian process spec and data spec */

  logf.file_name = *ap++;

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  gp_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  gp = logg.data['P'];
  m  = logg.data['M'];
  sv = logg.data['V'];

  gp_check_specs_present(gp,0,m,sv);

  h->total_hypers = gp_hyper_count(gp,m);

  logg.req_size['S'] = h->total_hypers * sizeof(double);

  if ((data_spec = logg.data['D'])==0)
  { fprintf(stderr,"No data specification in log file\n");
    exit(1);
  }

  M_targets = m!=0 && m->type=='C' ? gp->N_outputs : data_spec->N_targets;

  if (test_inputs_spec[0]!=0)  
  { strcpy (data_spec->test_inputs,  test_inputs_spec);
    strcpy (data_spec->test_targets, test_targets_spec);
  }

  have_targets = data_spec->test_targets[0]!=0;

  if (data_spec->test_inputs[0]==0)
  { fprintf(stderr,"No test inputs specified\n");
    exit(1);
  }

  /* Can't handle survival models. */

  if (m!=0 && m->type=='V')
  { fprintf(stderr,"Can't handle survival models in gp-pred\n");
    exit(1);
  }

  /* Check for illegal option combinations. */

  if (op_r && (m!=0 && m->type=='B' || m!=0 && m->type=='C')
   || op_p && (m==0) || op_l && (m==0 || m->type!='R')
   || op_m && (m==0 || m->type=='R' || m->type=='V')
   || op_d && (m!=0 && m->type=='B' || m!=0 && m->type=='C'))
  { fprintf(stderr,"Illegal combination of options with data model\n");
    exit(1);
  }

  if (op_b && op_a)
  { fprintf(stderr,"Options 'a' and 'b' are incompatible\n");
    exit(1);
  }

  if (!have_targets && (op_t || op_p || op_a))
  { fprintf(stderr,
      "Options 't', 'p', and 'a' are valid only when the targets are known\n");
    exit(1);
  }

  /* Read training & test data. */

  gp_data_read (1, 1, gp, m, sv);

  /* Allocate some space. */

  train_cov = chk_alloc (N_train*N_train, sizeof(double));

  scr1 = chk_alloc (N_train, sizeof(double));
  scr2 = chk_alloc (N_train, sizeof(double));
  scr3 = chk_alloc (N_train, sizeof(double));

  meanp = chk_alloc (N_test*gp->N_outputs, sizeof(double));
  varp  = chk_alloc (N_test*gp->N_outputs, sizeof(double));

  /* Look at things in order to make guesses and/or probability judgements. */

  if (guessing)
  {
    /* Allocate space for error accumulators. */
  
    log_prob     = chk_alloc (N_test, sizeof (double));
    sum_targets  = chk_alloc (M_targets*N_test, sizeof (double));
    prb_dist1    = chk_alloc (M_targets, sizeof (double));
    prb_dist2    = chk_alloc (M_targets, sizeof (double));
    sq_error     = chk_alloc (M_targets, sizeof (double));
    sq_error_sq  = chk_alloc (M_targets, sizeof (double));
    guess        = chk_alloc (M_targets, sizeof (double));
    error        = chk_alloc (M_targets, sizeof (double));
    abs_error    = chk_alloc (data_spec->N_targets, sizeof (double));
    wrong        = chk_alloc (data_spec->N_targets, sizeof (double));
    abs_error_sq = chk_alloc (data_spec->N_targets, sizeof (double));
    wrong_sq     = chk_alloc (data_spec->N_targets, sizeof (double));

    if (op_d)
    { 
      median_sample = chk_alloc (N_test, sizeof *median_sample);
  
      for (i = 0; i<N_test; i++)
      { 
        median_sample[i] = chk_alloc (data_spec->N_targets, 
                                      sizeof **median_sample);

        for (j = 0; j<data_spec->N_targets; j++)
        { median_sample[i][j] = chk_alloc (Median_sample*Max_median_points,
                                           sizeof ***median_sample);
        }
      }
    }

    /* Initialize various accumulators. */
  
    for (i = 0; i<N_test; i++)
    { for (j = 0; j<M_targets; j++)
      { sum_targets[M_targets*i+j] = 0;
      }
    }

    ms_count = 0;
  
    /* Go through all the iterations, taken from all the log files. */
  
    N_gps = 0;
  
    for (;;)
    {
      /* Figure out range to use from this log file. */

      if (**ap == '@') /* Range gives indexes via cpu times */
      {
        parse_time_range ((*ap++) + 1, &dlindex, &dhindex, &mod_no);

        lindex = -1;

        while (!logf.at_end && lindex==-1)
        { 
          log_gobble(&logf,&logg);

          if (logg.last_index>=0 && logg.data['i']!=0
           && logg.index['i']==logg.last_index)
          { it = logg.data['i'];
            if (it->time>=60000*dlindex)
            { lindex = logg.last_index;
            }
          }
        }  

        if (lindex==-1)
        { fprintf(stderr,"Warning: No records in range in %s\n",logf.file_name);
          goto close_file;
        }

        if (dhindex==-2) 
        { hindex = lindex;
        }
        else if (dhindex==-1)
        { log_file_last(&logf); 
          hindex = logf.header.index; 
        }
        else 
        {
          for (;;)
          { 
            if (logg.index['i']==logg.last_index)
            { it = logg.data['i'];
              if (it->time>60000*dhindex)
              { break;
              }
              hindex = logg.last_index;
            }

            if (logf.at_end) break;

            log_gobble(&logf,&logg);
          }
        }

        log_file_first(&logf);
      }

      else /* Range doesn't start with "@", gives indexes directly */
      {
        parse_range (*ap++, &lindex, &hindex, &mod_no);

        if (hindex==-2) 
        { hindex = lindex;
        }
        else if (hindex==-1) 
        { log_file_last(&logf);
          hindex = logf.header.index;
          log_file_first(&logf); 
        }
      }

      /* Reduce number of iterations asked for if it's more than is available.*/

      if (mod_no<0 && hindex-lindex+1<-mod_no)
      { mod_no = - (hindex-lindex+1);
      }

      /* Read iterations in this range. */

      if (0) /* Debug message */
      { fprintf (stderr, "Using index range of %d to %d for log file %s.\n", 
          lindex, hindex, logf.file_name);
      }

      if (mod_no<0) target_index = lindex;

      for (;;)
      {
        /* Skip to next desired index, or to end of range. */

        while (!logf.at_end && logf.header.index<=hindex
         && (logf.header.index<lindex || 
             (mod_no>0 ? logf.header.index%mod_no!=0
                       : logf.header.index<(int)(target_index+0.5))))
        { 
          log_file_forward(&logf);
        }
  
        if (logf.at_end || logf.header.index>hindex)
        { break;
        }
    
        /* Gobble up records for this index. */
  
        log_gobble(&logf,&logg);
       
        /* Apply the Gaussian process with this index to the test cases. */

        if (0) /* Debug message */
        { fprintf (stderr, "Using record from log file %s at index %d.\n",
            logf.file_name, logg.last_index);
        }
    
        if (logg.index['S']!=logg.last_index)
        { fprintf(stderr,
                  "Warning: Missing data for iteration with index %d in %s\n",
                  logg.last_index, logf.file_name);
        }
        else
        { 
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

          if (op_d && N_gps>=Max_median_points)
          { fprintf(stderr,
              "Too many iterations being used to find median (max %d)\n",
              Max_median_points);
            exit(1);
          }

          /* Find predictive means and (maybe) variances for each output
             of every test case, storing them in meanp and varp. */

          for (j = 0; j<gp->N_outputs; j++)
          {
            /* Compute covariance matrix for training cases, then invert.
               Don't redo computations if previous matrix is still applicable.*/

            if (j==0 || m!=0 && m->type=='R' && !have_values
              && (m->noise.alpha[2]!=0 || *hypers.noise[j]!=*hypers.noise[j-1]))
            {
              gp_train_cov (gp, m, h, j, noise_variances, 
                            have_values ? train_cov : 0,
                            have_values ? 0 : train_cov,
                            0);
        
              if (!cholesky (train_cov, N_train, 0)
               || !inverse_from_cholesky (train_cov, scr1, scr2, N_train))
              { fprintf (stderr,
                 "Couldn't invert covariance matrix of training cases\n");
                exit(1);
              }
            }

            /* Pre-compute product of inverse covariance for training points
               with training targets or latent values, for use in finding 
               predictive means. */

            for (i = 0; i<N_train; i++)
            { scr1[i] = inner_product (train_cov+i*N_train, 1, 
                           have_values ? latent_values+j : train_targets+j, 
                           gp->N_outputs, N_train);
            }

            /* Find the predictive mean/variance for this output for each 
               test case. */

            for (i = 0; i<N_test; i++)
            {
              /* Find covariances between training cases and test case. */

              gp_cov (gp, h, test_inputs+gp->N_inputs*i, 1, 
                      train_inputs, N_train, scr2, 0);

              /* Find predictive mean using the vector precomputed above. */

              meanp[i*gp->N_outputs+j] 
                = inner_product (scr1, 1, scr2, 1, N_train);

              /* Find the predictive variance if it will be needed. */

              if (op_p || op_d || op_r && data_spec->target_trans[j].take_log
               || m!=0 && m->type!='R')
              { 
                gp_cov (gp, h, test_inputs+gp->N_inputs*i, 1,
                               test_inputs+gp->N_inputs*i, 1, 
                               &varp[i*gp->N_outputs+j], 0);

                if (gp->has_jitter)
                { varp[i*gp->N_outputs+j] += exp(2 * *h->jitter);
                }

                matrix_product (scr2, train_cov, scr3, 1, N_train, N_train);

                varp[i*gp->N_outputs+j] 
                  -= inner_product (scr3, 1, scr2, 1, N_train);

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
            rand_seed(101*logg.last_index+i);

            if (m==0 || m->type=='R') /* Model for real data */
            {
              curr_log_prob = 0;

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
                  tr = &data_spec->target_trans[j];
  
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

                sum_targets[M_targets*i+j] += meanp[i*gp->N_outputs+j];
              
                if (have_targets && op_p) 
                { 
                  targ = test_targets[gp->N_outputs*i+j];

                  curr_log_prob += - log(2*M_PI*varp[i*gp->N_outputs+j])/2 
                    - ((targ-meanp[i*gp->N_outputs+j])
                       * (targ-meanp[i*gp->N_outputs+j])) 
                         / (2*varp[i*gp->N_outputs+j]);
  
                  if (op_r)
                  { tr = &data_spec->target_trans[j];
                    curr_log_prob += log(tr->scale);
                    if (tr->take_log)
                    { curr_log_prob -= log (data_inv_trans 
                        (test_targets[data_spec->N_targets*i+j], *tr));
                    }
                  }
                }
  
                if (op_d)
                {
                  for (k = 0; k<Median_sample; k++)
                  {
                    guessp = meanp[i*gp->N_outputs+j] 
                              + sqrt(varp[i*gp->N_outputs+j])*rand_gaussian();
  
                    if (op_r)
                    { guessp = data_inv_trans (guessp, 
                                               data_spec->target_trans[j]);
                    }

                    median_sample[i][j][ms_count+k] = guessp;
                  }
                }
              }
            }

            else if (m->type=='B') /* Model for binary data */
            {
              curr_log_prob = 0;
 
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

                sum_targets[M_targets*i+j] += prb;

                if (have_targets && op_p)
                { targ = test_targets[gp->N_outputs*i+j];
                  curr_log_prob += targ==1 ? log(prb) : log(1-prb);
                }
              }
            }

            else if (m->type=='C') /* Model for multi-class data */
            {
              curr_log_prob = 0;
 
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
              { sum_targets[M_targets*i+j] += prb_dist1[j];
              }

              if (have_targets && op_p)
              { targ = test_targets[i];
                curr_log_prob += log(prb_dist1[(int)targ]);
              }
            }

            else /* Type of model that isn't handled */
            { 
              abort();
            }

            if (have_targets && op_p)
            {
              log_prob[i] = N_gps==0 ? curr_log_prob
                                     : addlogs(log_prob[i],curr_log_prob);
            }
          }

          ms_count += Median_sample;          
          N_gps += 1;
        }


        /* Find next target index, if going by number of networks. */

        if (mod_no<0) target_index += (double) (hindex-lindex) / (-mod_no-1);
      }

    close_file:  
      log_file_close(&logf);
  
      /* See if we're done. */
  
      if (*ap==0) break;
  
      /* Open next log file.  Note that initial records (architecture, etc.) 
         are read only from the first log file. */
  
      logf.file_name = *ap++;
      log_file_open(&logf,0);
  
      logg.index['S'] = -1;  /* So they won't be mistaken for records from */
      logg.index['V'] = -1;  /* the new log file.                          */
      logg.index['F'] = -1;
    }
  
    if (N_gps==0) 
    { fprintf(stderr,"None of the specified iterations were found\n");
      exit(1);
    }
  }

  /* Print headings. */

  if (!op_b) printf("\nNumber of iterations used: %d\n",N_gps);

  if (!op_b && !op_a)
  { 
    printf("\nCase");

    if (op_i)
    { printf(" %-*s",7*gp->N_inputs," Inputs");
    }

    if (op_t)
    { printf(" %-*s",7*data_spec->N_targets," Targets");
    }

    if (op_p)
    { printf("  Log Prob");
    }

    if (op_m)
    { printf(" %-*s",3*data_spec->N_targets," Guesses");
      if (have_targets) 
      { printf(" %-*s",2*data_spec->N_targets," Wrong?");
      }
    }

    if (op_n)
    { printf(" %-*s",7*M_targets," Means");
      if (have_targets) 
      { printf(" %-*s",7*data_spec->N_targets," Error^2");
      }
    }

    if (op_d)
    { printf(" %-*s",7*M_targets," Medians");
      if (have_targets) 
      { printf(" %-*s",7*M_targets," |Error|");
      }
    }

    printf("\n\n");
  }

  /* Make any required guesses and print the results for each case. */

  if (guessing)
  {
    ave_log_prob = ave_log_prob_sq = 0;
 
    for (j = 0; j<M_targets; j++) 
    { sq_error[j] = 0;
      sq_error_sq[j] = 0;
    }
 
    for (j = 0; j<data_spec->N_targets; j++) 
    { abs_error[j] = wrong[j] = 0;
      abs_error_sq[j] = wrong_sq[j] = 0;
    }

    tsq_error = tabs_error =  0;
    tsq_error_sq = tabs_error_sq = 0;
  }

  for (i = 0; i<N_test; i++)
  { 
    if (put_in_blanks)
    { if (i>0 
       && test_inputs[gp->N_inputs*i+0] != test_inputs[gp->N_inputs*(i-1)+0])
      { printf("\n");
      }
    }

    if (!op_a && !op_b) printf("%4d",i+1);
 
    if (!op_a && op_i)
    { printf(" ");
      for (j = m!=0 && m->type=='V' && sv->hazard_type!='C' ? 1 : 0; 
               j<gp->N_inputs; j++)
      { printf(op_b ? " %+.8e" : " %6.2f",test_inputs[gp->N_inputs*i+j]);
      }
    }

    if (!op_a && op_t)
    { printf(" ");
      for (j = 0; j<data_spec->N_targets; j++)
      { double val;
        val = test_targets[data_spec->N_targets*i+j];
        if (op_r) val = data_inv_trans(val,data_spec->target_trans[j]);
        printf(op_b ? " %+.8e" : " %6.2f",val);
      }
      if (data_spec->N_targets==1) printf(" ");
    }

    if (op_p)
    { log_prob[i] -= log((double)N_gps);
      ave_log_prob += log_prob[i];
      ave_log_prob_sq += log_prob[i]*log_prob[i];
      if (!op_a) printf(op_b ? " %.8e" : " %9.3f",log_prob[i]);
    }

    if (op_m)
    { 
      if (m->type=='C')
      { int g;
        g = 0;
        for (j = 1; j<gp->N_outputs; j++)
        { if (sum_targets[M_targets*i+j] > sum_targets[M_targets*i+g])
          { g = j;
          }
        }
        guess[0] = g;
        if (have_targets)
        { error[0] = guess[0] != test_targets[data_spec->N_targets*i];
          wrong[0] += error[0];
          wrong_sq[0] += error[0]*error[0];
        }
      }

      if (m->type=='B')
      { for (j = 0; j<gp->N_outputs; j++)
        { guess[j] = sum_targets[M_targets*i+j]/N_gps > 0.5;
          if (have_targets)
          { error[j] = guess[j] != test_targets[data_spec->N_targets*i+j];
            wrong[j] += error[j];
            wrong_sq[j] += error[j]*error[j];
          }
        }
      }

      if (!op_a) 
      { printf(" ");
        if (data_spec->N_targets<2) printf("   ");
        for (j = 0; j<data_spec->N_targets; j++) 
        { printf(" %2.0f",guess[j]);
        }
        if (data_spec->N_targets<3) printf("  ");
        if (have_targets)
        { printf(" ");
          if (data_spec->N_targets<2) printf("   ");
          for (j = 0; j<data_spec->N_targets; j++) 
          { printf(" %1.0f",error[j]);
          }
        }
      }
    }

    if (op_n)
    { 
      double tot_error;

      tot_error = 0;

      for (j = 0; j<M_targets; j++)
      { guess[j] = sum_targets[M_targets*i+j] / N_gps;
        if (have_targets)
        { double val, diff;
          val = m!=0 && m->type=='C' ? j==test_targets[i]
                                     : test_targets[data_spec->N_targets*i+j];
          if (op_r) val = data_inv_trans(val,data_spec->target_trans[j]);
          if (m!=0 && m->type=='V' && val<0) val = -val;
          diff = guess[j] - val;
          error[j] = diff * diff;
          tot_error += error[j];
          sq_error[j] += error[j];
          sq_error_sq[j] += error[j]*error[j];
        }
      }

      tsq_error += tot_error;
      tsq_error_sq += tot_error * tot_error;

      if (!op_a) 
      { printf(" ");
        for (j = 0; j<M_targets; j++) 
        { printf(op_b ? " %+.8e" : " %6.2f",guess[j]);
        }
        if (have_targets)
        { printf(" ");
          if (m!=0 && m->type=='C') 
          { printf(op_b ? " %.8e" : " %6.4f",tot_error);
          }
          else
          { for (j = 0; j<data_spec->N_targets; j++) 
            { printf(op_b ? " %.8e" : " %6.4f",error[j]);
            }
          }
        }
      }
    }

    if (op_d)
    { 
      double tot_error;

      tot_error = 0;

      for (j = 0; j<data_spec->N_targets; j++)
      { guess[j] = find_median (median_sample[i][j], ms_count);
        if (have_targets)
        { double val, diff;
          val = test_targets[data_spec->N_targets*i+j];
          if (op_r) val = data_inv_trans(val,data_spec->target_trans[j]);
          if (m!=0 && m->type=='V' && val<0) val = -val;
          diff = guess[j] - val;
          error[j] = diff<0 ? -diff : diff;
          tot_error += error[j];
          abs_error[j] += error[j];
          abs_error_sq[j] += error[j]*error[j];
        }
      }

      tabs_error += tot_error;
      tabs_error_sq += tot_error * tot_error;
  
      if (!op_a)
      { printf(" ");
        if (data_spec->N_targets==1) printf(" ");
        for (j = 0; j<data_spec->N_targets; j++) 
        { printf(op_b ? " %+.8e" : " %6.2f",guess[j]);
        }
        if (have_targets)
        { printf(" ");
          for (j = 0; j<data_spec->N_targets; j++) 
          { printf(op_b ? " %.8e" : " %6.4f",error[j]);
          }
        }
      }
    }
 
    if (!op_a) printf("\n");
  }

  /* Print the averages. */

  if (have_targets && guessing && !op_b)
  { 
    printf("\n");

    if (op_a)
    { printf("Number of test cases: %d\n\n", N_test);
    }

    if (op_p) 
    { double a;
      a = ave_log_prob/N_test;
      printf("Average log probability of targets: %9.3f", a);
      if (N_test>1) 
      { printf("+-%.3f", sqrt((ave_log_prob_sq/N_test-a*a) / (N_test-1)));
      }
      printf("\n\n");
    }

    if (op_m) 
    { double a;
      printf("Fraction of guesses that were wrong: ");
      for (j = 0; j<data_spec->N_targets; j++)
      { a = wrong[j] / N_test;
        printf(" %6.4f",a);
        if (N_test>1) 
        { printf("+-%6.4f", sqrt((wrong_sq[j]/N_test-a*a) / (N_test-1)));
        }
      }
      printf("\n");
    }

    if (op_n) 
    { 
      double a;

      printf("Average squared error guessing mean: ");

      if (m!=0 && m->type=='C')
      { a = tsq_error / N_test;
        printf (" %8.5f", a);
        if (N_test>1) 
        { printf("+-%.5f",sqrt((tsq_error_sq/N_test - a*a) / (N_test-1)));
        }
        printf("\n");
      } 
      else
      { for (j = 0; j<data_spec->N_targets; j++)
        { a = sq_error[j]/N_test;
          printf (" %8.5f", a);
          if (N_test>1)
          { printf("+-%.5f", sqrt((sq_error_sq[j]/N_test - a*a) / (N_test-1)));
          }
        }
        printf("\n");
      }

      if (data_spec->N_targets>1) 
      { a = tsq_error / N_test;
        printf("                                            (total %.5f",a);
        if (N_test>1) 
        { printf("+-%.5f",sqrt((tsq_error_sq/N_test - a*a) / (N_test-1)));
        }
        printf(")\n");
      }
    }

    if (op_d) 
    { 
      double a;

      printf("Average abs. error guessing median:  ");

      for (j = 0; j<data_spec->N_targets; j++)
      { a = abs_error[j] / N_test;
        printf(" %8.5f", a);
        if (N_test>1)
        { printf("+-%.5f", sqrt((abs_error_sq[j]/N_test - a*a) / (N_test-1)));
        }
      }
      printf("\n");

      if (data_spec->N_targets>1)
      { a = tabs_error / N_test;
        printf("                                            (total %.5f",a);
        if (N_test>1) 
        { printf("+-%.5f",sqrt((tabs_error_sq/N_test - a*a) / (N_test-1)));
        }
        printf(")\n");
      }
    }
  }

  if (!op_b)    
  { printf("\n");
  }
  
  exit(0);
}


/* FIND MEDIAN OF ARRAY OF VALUES. */

static int flt_cmp (const void *a, const void *b) 
{ return *(float*)a > *(float*)b ? +1 : *(float*)a < *(float*)b ? -1 : 0;
}

static float find_median
( float *data,		/* Array of data values */
  int n			/* Number of values in array */
)
{
  qsort (data, n, sizeof *data, flt_cmp);
 
  return n%2==1 ? data[(int)(n/2)] 
                : (data[(int)(n/2)] + data[(int)((n-1)/2)]) / 2;
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf(stderr,
   "Usage: gp-pred options { log-file range } [ / test-inputs [ test-targets ] ]\n");
  fprintf(stderr,
    "  Opt: i=inputs t=targets p=prob m=mode n=mean d=median\n");
  fprintf(stderr,
    "       r = raw b/B=bare a=aveonly\n");
  exit(1);
}
