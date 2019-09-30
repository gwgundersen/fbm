/* PRED.C - Skeleton for programs that make predictions for test cases. */

/* Copyright (c) 1995, 1996, 1998, 1999 by Radford M. Neal 
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
#include "log.h"
#include "prior.h"
#include "model.h"
#include "data.h"
#include "mc.h"
#include "pred.h"


/* DEBUG FLAG.  Normally set to 0. */

#define DEBUG 0


/* SHARED VARIABLES.  Declared as external in pred.h, which is where they're
   documented. */

int op_i, op_t, op_r, 
    op_p, op_m, op_n, 
    op_d, op_l, op_b, 
    op_a, op_q, op_N;

int keep[10];

log_gobbled logg;

model_specification *m;
model_survival *sv;

data_transformation *tr;

int M_targets;

int N_records_used;
double *test_targ_pred;
double *test_log_prob;

float ***median_sample;
int ms_count;
 
int have_targets;

int use_inverse;
int alt_mean;


/* LOCAL PROCEDURES. */

static float find_median (float *, int);

static void usage(void);


/* MAIN PROGRAM. */

void main
( int argc,
  char **argv
)
{
  double *sq_error, *sq_error_sq, *abs_error, *abs_error_sq, *wrong, *wrong_sq;
  double tsq_error, tsq_error_sq, tabs_error, tabs_error_sq;
  double ave_log_prob, ave_log_prob_sq;
  double *guess, *error;
  double *sum_targets;
  double *log_prob;

  char test_inputs_spec[100];
  char test_targets_spec[100];

  log_file logfile;

  double dlindex, dhindex, target_index;
  int lindex, hindex, mod_no;

  int N_records_used, has_weights;
  double max_log_weight, sum_weights;
  double sumsq_weights;

  mc_temp_state *ts;
  mc_iter *it;
  char *options;
  int put_in_blanks;
  int guessing;

  double lf, f, w;
  int i, j, l;
  char **ap;

  /* Special fudge to handle computational options for gp-pred. */

  if (strcmp(pred_app_name,"gp")==0)
  { for (;;)
    { if (argc>1 && strcmp(argv[1],"-use-inverse")==0)
      { use_inverse = 1;
      }
      else if (argc>1 && strcmp(argv[1],"-alt-mean")==0)
      { alt_mean = 1;
      }
      else
      { break;
      }
      argc -= 1;
      argv += 1;
    }
  }

  if (alt_mean && use_inverse)
  { fprintf(stderr,
      "Can't use -alt-mean and -use-inverse options at the same time\n");
    exit(1);
  }

  /* Look at other arguments. */

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
  op_q = strchr(options,'q')!=0;
  op_l = strchr(options,'l')!=0;
  op_b = strchr(options,'b')!=0 || strchr(options,'B')!=0;
  op_a = strchr(options,'a')!=0;


  op_N = 0;
  for (l = 0; l<10; l++) keep[l] = 0;
  for (i = 0; options[i]!=0; i++)
  { if (options[i]>='0' && options[i]<='9')
    { op_N += 1;
      if (keep[options[i]-'0']) usage();
      keep[options[i]-'0'] = 1;
    }
  }

  if (strlen(options) 
       != op_i+op_t+op_r+op_p+op_m+op_n+op_d+op_b+op_a+op_q+op_l+op_N)
  { usage();
  }

  guessing = op_p || op_m || op_n || op_d;

  put_in_blanks = strchr(options,'B')!=0;

  if ((argc-2)%2!=0 && (guessing || argc!=3)) usage();

  ap = argv+2;

  /* Open first log file and read specifications. */

  logfile.file_name = *ap++;

  log_file_open (&logfile, 0);

  log_gobble_init(&logg,0);
  pred_app_record_sizes();

  while (!logfile.at_end && logfile.header.index<0)
  { log_gobble(&logfile,&logg);
  }

  m  = logg.data['M'];
  sv = logg.data['V'];

  if ((data_spec = logg.data['D'])==0)
  { fprintf(stderr,"No data specification in log file\n");
    exit(1);
  }

  M_targets = m!=0 && m->type=='C' ? data_spec->int_target 
                                   : data_spec->N_targets;

  if (test_inputs_spec[0]!=0)  
  { strcpy (data_spec->test_inputs,  test_inputs_spec);
    strcpy (data_spec->test_targets, test_targets_spec);
  }

  have_targets = data_spec->test_targets[0]!=0;

  if (data_spec->test_inputs[0]==0)
  { fprintf(stderr,"No test inputs specified\n");
    exit(1);
  }

  /* Check for illegal option combinations. */

  if (op_a && (op_t || op_i || op_q)
   || op_r && (m!=0 && m->type=='B' || m!=0 && m->type=='C')
   || op_p && (m==0) || op_l && (m==0 || m->type!='R')
   || op_m && (m==0 || m->type=='R' || m->type=='V')
   || (op_d || op_q) && (m!=0 && m->type=='B' || m!=0 && m->type=='C'))
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

  /* Initialize application, including reading data. */

  pred_app_init();

  /* Look at things in order to make guesses and/or probability judgements. */

  if (argc!=3)
  {
    /* Allocate space for error accumulators. */
  
    log_prob       = chk_alloc (N_test, sizeof (double));
    sum_targets    = chk_alloc (M_targets*N_test, sizeof (double));
    test_log_prob  = chk_alloc (N_test, sizeof (double));
    test_targ_pred = chk_alloc (M_targets*N_test, sizeof (double));
    sq_error       = chk_alloc (M_targets, sizeof (double));
    sq_error_sq    = chk_alloc (M_targets, sizeof (double));
    guess          = chk_alloc (M_targets, sizeof (double));
    error          = chk_alloc (M_targets, sizeof (double));
    abs_error      = chk_alloc (data_spec->N_targets, sizeof (double));
    wrong          = chk_alloc (data_spec->N_targets, sizeof (double));
    abs_error_sq   = chk_alloc (data_spec->N_targets, sizeof (double));
    wrong_sq       = chk_alloc (data_spec->N_targets, sizeof (double));

    if (op_d || op_q)
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
  
    N_records_used = 0;
    has_weights = 0;
    sum_weights = 0;
    sumsq_weights = 0;
  
    for (;;)
    {
      /* Figure out range to use from this log file. */

      if (**ap == '@') /* Range gives indexes via cpu times */
      {
        parse_time_range ((*ap++) + 1, &dlindex, &dhindex, &mod_no);

        lindex = -1;

        while (!logfile.at_end && lindex==-1)
        { 
          log_gobble(&logfile,&logg);

          if (logg.last_index>=0 && logg.data['i']!=0
           && logg.index['i']==logg.last_index)
          { it = logg.data['i'];
            if (it->time>=60000*dlindex)
            { lindex = logg.last_index;
            }
          }
        }  

        if (lindex==-1)
        { fprintf(stderr,
                  "Warning: No records in range in %s\n",logfile.file_name);
          goto close_file;
        }

        if (dhindex==-2) 
        { hindex = lindex;
        }
        else if (dhindex==-1)
        { log_file_last(&logfile); 
          hindex = logfile.header.index; 
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

            if (logfile.at_end) break;

            log_gobble(&logfile,&logg);
          }
        }

        log_file_first(&logfile);
      }

      else /* Range doesn't start with "@", gives indexes directly */
      {
        parse_range (*ap++, &lindex, &hindex, &mod_no);

        if (hindex==-2) 
        { hindex = lindex;
        }
        else if (hindex==-1) 
        { log_file_last(&logfile);
          hindex = logfile.header.index;
          log_file_first(&logfile); 
        }
      }

      /* Reduce number of iterations asked for if it's more than is available.*/

      if (mod_no<0 && hindex-lindex+1<-mod_no)
      { mod_no = - (hindex-lindex+1);
      }

      /* Read iterations in this range. */

      if (DEBUG)
      { fprintf (stderr, "Using index range of %d to %d for log file %s.\n", 
          lindex, hindex, logfile.file_name);
      }

      if (mod_no<0) target_index = lindex;

      for (;;)
      {
        /* Skip to next desired index, or to end of range. */

        while (!logfile.at_end && logfile.header.index<=hindex
         && (logfile.header.index<lindex || 
             (mod_no>0 ? logfile.header.index%mod_no!=0
                       : logfile.header.index<(int)(target_index+0.5))))
        { 
          log_file_forward(&logfile);
        }
  
        if (logfile.at_end || logfile.header.index>hindex)
        { break;
        }
    
        /* Gobble up records for this index. */
  
        log_gobble(&logfile,&logg);

        if (DEBUG)
        { fprintf (stderr, "Using record from log file %s at index %d.\n",
            logfile.file_name, logg.last_index);
        }

        if ((op_d || op_q) && N_records_used>=Max_median_points)
        { fprintf(stderr,
           "Too many iterations being used to find median/quantiles (max %d)\n",
            Max_median_points);
          exit(1);
        }

        /* Look at weight in iteration record and tempering state. */

        it = logg.data['i'] != 0 && logg.index['i']==logg.last_index
               ? logg.data['i'] : 0;
        ts  = logg.data['b']!=0 && logg.index['b']==logg.last_index
               ? logg.data['b'] : 0;

        if (it!=0 && it->log_weight!=0)
        { if (op_d || op_q)
          { fprintf(stderr,
           "Median/quantiles are not implemented for weighted data from AIS\n");
            exit(1);
          }
          has_weights = 1;
          w = it->log_weight;
        }
        else
        { w = 0;
        }
       
        /* Use records at this index to make predictions for the test cases. */

        if (ts!=0 && ts->inv_temp!=1)
        { /* Ignore - wrong temperature */
        }
        else if (!pred_app_use_index())
        { fprintf(stderr,
            "Warning: Missing data at index %d in %s - ignored this index\n",
             logg.last_index, logfile.file_name);
        }
        else
        { 
          if (N_records_used==0)
          { max_log_weight = w;
            lf = 0;
            f = 1;
          }
          else if (w>max_log_weight)
          { lf = max_log_weight-w;
            f = exp(lf);
            sum_weights *= f;
            sumsq_weights *= f*f;
            if (have_targets)
            { for (i = 0; i<N_test; i++) 
              { log_prob[i] += lf;
              }
            }
            for (j = N_test*M_targets-1; j>=0; j--) 
            { sum_targets[j] *= f;
            }
            max_log_weight = w;
            lf = 0;
            f = 1;
          }
          else
          { lf = w-max_log_weight;
            f = exp(lf);
          }

          sum_weights += f;
          sumsq_weights += f*f;

          if (have_targets)
          { for (i = 0; i<N_test; i++)
            { log_prob[i] = N_records_used==0 ? test_log_prob[i] 
                : addlogs (log_prob[i], test_log_prob[i]+lf);
            }
          }

          for (j = N_test*M_targets-1; j>=0; j--) 
          { sum_targets[j] += f * test_targ_pred[j];
          }

          ms_count += Median_sample; 
          N_records_used += 1;
        }

        /* Find next target index, if going by number of records. */

        if (mod_no<0) target_index += (double) (hindex-lindex) / (-mod_no-1);
      }

    close_file:  
      log_file_close(&logfile);
  
      /* See if we're done. */
  
      if (*ap==0) break;
  
      /* Open next log file.  Note that initial records (architecture, etc.) 
         are read only from the first log file. */
  
      logfile.file_name = *ap++;
      log_file_open(&logfile,0);

      pred_app_finish_file();  
    }
  
    if (N_records_used==0) 
    { fprintf(stderr,"None of the specified iterations were found\n");
      exit(1);
    }
  }

  /* Print overall information. */

  if (!op_b && argc!=3)
  { printf("\nNumber of iterations used: %d\n",N_records_used);
    if (has_weights && N_records_used>1)
    { printf("Adjusted sample size: %.1f\n",
        sum_weights*sum_weights / sumsq_weights);
    }
    if (has_weights)
    { printf("\nLog of marginal likelihood: %.3f",
        log(sum_weights/N_records_used) + max_log_weight);
      if (N_records_used>1 && sum_weights*sum_weights/sumsq_weights>=2)
      { printf(" +- %.3f",
          sqrt (sumsq_weights/(sum_weights*sum_weights) - 1.0/N_records_used));
      }
      else
      { printf(" +- ???");
      }
      printf("\n");
    }
  }

  /* Print headings. */

  if (!op_b && !op_a)
  { 
    printf("\nCase");

    if (op_i)
    { printf(" %-*s",7*data_spec->N_inputs," Inputs");
    }

    if (op_t)
    { printf(" %-*s",7*data_spec->N_targets,"Targets");
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
    { printf(" %-*s",7*M_targets,"  Means");
      if (have_targets) 
      { printf(" %-*s",7*data_spec->N_targets,"Error^2");
      }
    }

    if (op_d)
    { printf(" %-*s",7*M_targets,"Medians");
      if (have_targets) 
      { printf(" %-*s",7*M_targets,"|Error|");
      }
    }

    if (op_q)
    { printf(" %-*s",7*M_targets,"10% Qnt");
      printf(" %-*s",7*M_targets,"90% Qnt");
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
    { if (i>0 && test_inputs[data_spec->N_inputs*i+0] 
                  != test_inputs[data_spec->N_inputs*(i-1)+0])
      { printf("\n");
      }
    }

    if (!op_a && !op_b) printf("%4d",i+1);
 
    if (!op_a && op_i)
    { printf(" ");
      for (j = 0; j<data_spec->N_inputs; j++)
      { printf(op_b ? " %+.8e" : " %6.2f",test_inputs[data_spec->N_inputs*i+j]);
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
    }

    if (op_p)
    { log_prob[i] -= log(sum_weights);
      ave_log_prob += log_prob[i];
      ave_log_prob_sq += log_prob[i]*log_prob[i];
      if (!op_a) printf(op_b ? " %.8e" : " %9.3f",log_prob[i]);
    }

    if (op_m)
    { 
      if (m->type=='C')
      { int g;
        g = 0;
        for (j = 1; j<M_targets; j++)
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
      { for (j = 0; j<data_spec->N_targets; j++)
        { guess[j] = sum_targets[M_targets*i+j]/sum_weights > 0.5;
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
      { guess[j] = sum_targets[M_targets*i+j] / sum_weights;
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

    if (op_q)
    { 
      for (j = 0; j<data_spec->N_targets; j++)
      { (void) find_median (median_sample[i][j], ms_count);
      }

      printf(" ");
      for (j = 0; j<data_spec->N_targets; j++) 
      { printf(op_b ? " %+.8e" : " %6.2f",
               median_sample[i][j][(int)(0.5 + 0.1*ms_count)]);
      }

      printf(" ");
      for (j = 0; j<data_spec->N_targets; j++) 
      { printf(op_b ? " %+.8e" : " %6.2f",
               median_sample[i][j][(int)(0.5 + 0.9*ms_count)]);
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


/* FIND MEDIAN OF ARRAY OF VALUES.  The array is sorted in increasing order
   in the process, which allows other quantiles to be found by looking at
   the appropriate element. */

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
   "Usage: %s-pred options { log-file range } [ / test-inputs [ test-targets ] ]\n", pred_app_name);
  fprintf(stderr,
    "  Opt: i=inputs, t=targets, p=prob, m=mode, n=mean, d=median, q=quantiles\n");
  fprintf(stderr,
    "       <digit> = select component, r = raw data, b/B = bare, a = averages only\n");
  exit(1);
}
