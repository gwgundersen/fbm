/* NET-PRED.C - Program to make predictions for test cases. */

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
 *
 * Modifications to allow selection of a given number of networks and to
 * allow specification of ranges by cpu-time are adapted from modifications
 * by Carl Edward Rasmussen, 1995 
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "misc.h"
#include "net.h"
#include "log.h"
#include "data.h"
#include "numin.h"
#include "net-data.h"
#include "mc.h"


#define Max_median_nets 200	/* Max networks that can be used for median */
#define Median_sample 11	/* Size of sample per network */

static float find_median (float *, int);

static void usage(void);


/* MAIN PROGRAM. */

void main
( int argc,
  char **argv
)
{
  net_arch   *a;
  net_priors *p;

  net_sigmas sigmas, *s = &sigmas;
  net_params params, *w = &params;

  mc_iter    *it;

  int op_i, op_t, op_r, op_p, op_m, op_n, op_d, op_b, op_a, op_z, op_N;
  int keep[10];
  int put_in_blanks;
  int guessing;
  char *options;

  int M_targets, N_nets;

  double *sum_targets;
  double *curr_targets;

  double *sq_error, *sq_error_sq, *abs_error, *abs_error_sq, *wrong, *wrong_sq;
  double tsq_error, tsq_error_sq, tabs_error, tabs_error_sq;
  double *guess, *error;

  double *log_prob, ave_log_prob, ave_log_prob_sq, curr_log_prob;

  float ***median_sample;
  int ms_count;

  log_file logf;
  log_gobbled logg;

  char test_inputs_spec[100];
  char test_targets_spec[100];
 
  int have_targets;

  data_transformation *tr;

  double dlindex, dhindex, target_index;
  int lindex, hindex, mod_no;
  char **ap;
  int i, j, k, l;

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
  op_b = strchr(options,'b')!=0 || strchr(options,'B')!=0;
  op_a = strchr(options,'a')!=0;
  op_z = strchr(options,'z')!=0;

  op_N = 0;
  for (l = 0; l<10; l++) keep[l] = 0;
  for (i = 0; options[i]!=0; i++)
  { if (options[i]>='0' && options[i]<='9')
    { op_N += 1;
      if (keep[options[i]-'0']) usage();
      keep[options[i]-'0'] = 1;
    }
  }

  if (strlen(options) != op_i+op_t+op_r+op_p+op_m+op_n+op_d+op_b+op_a+op_z+op_N)
  { usage();
  }

  guessing = op_p || op_m || op_n || op_d;

  put_in_blanks = strchr(options,'B')!=0;

  if ((argc-2)%2!=0 && (guessing || argc!=3)) usage();

  ap = argv+2;

  /* Open first log file and read network architecture, priors, and data spec */

  logf.file_name = *ap++;

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  logg.req_size['A'] = sizeof *a;
  logg.req_size['P'] = sizeof *p;
  logg.req_size['D'] = sizeof *data_spec;

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }
  
  if ((a = logg.data['A'])==0)
  { fprintf(stderr,"No architecture specification in log file\n");
    exit(1);
  }

  s->total_sigmas = net_setup_sigma_count(a);
  w->total_params = net_setup_param_count(a);

  logg.req_size['S'] = s->total_sigmas * sizeof(net_sigma);
  logg.req_size['W'] = w->total_params * sizeof(net_param);

  if ((p = logg.data['P'])==0)
  { fprintf(stderr,"No prior specification in log file\n");
    exit(1);
  }

  if ((data_spec = logg.data['D'])==0)
  { fprintf(stderr,"No data specification in log file\n");
    exit(1);
  }

  M_targets = a->data_model=='C' ? a->N_outputs : data_spec->N_targets;

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

  if (op_r && (a->data_model=='B' || a->data_model=='C')
   || op_p && (a->data_model==0)
   || op_m && (a->data_model=='R' || a->data_model==0)
   || op_n && (a->data_model=='C')
   || op_d && (a->data_model=='B' || a->data_model=='C'))
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

  /* Read test data. */

  net_data_read (0, 1, a);

  /* Look at network to make guesses and/or probability judgements. */

  if (guessing)
  {
    /* Allocate space for error accumulators. */
  
    sum_targets  = chk_alloc (M_targets*N_test, sizeof (double));
    curr_targets = chk_alloc (M_targets, sizeof (double));
    log_prob     = chk_alloc (N_test, sizeof (double));
    guess        = chk_alloc (data_spec->N_targets, sizeof (double));
    error        = chk_alloc (data_spec->N_targets, sizeof (double));
    sq_error     = chk_alloc (data_spec->N_targets, sizeof (double));
    abs_error    = chk_alloc (data_spec->N_targets, sizeof (double));
    wrong        = chk_alloc (data_spec->N_targets, sizeof (double));
    sq_error_sq  = chk_alloc (data_spec->N_targets, sizeof (double));
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
        { median_sample[i][j] = chk_alloc (Median_sample*Max_median_nets,
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
  
    /* Go through all the networks, taken from all the log files. */
  
    N_nets = 0;
  
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

      /* Reduce number of networks asked for if it's more than is available. */

      if (mod_no<0 && hindex-lindex+1<-mod_no)
      { mod_no = - (hindex-lindex+1);
      }

      /* Read networks in this range. */

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
       
        /* Apply the network with this index to the test cases. */

        if (0) /* Debug message */
        { fprintf (stderr, "Using record from log file %s at index %d.\n",
            logf.file_name, logg.last_index);
        }
    
        if (logg.index['W']!=logg.last_index 
         || logg.index['S']!=logg.last_index)
        { fprintf(stderr,
                  "Warning: Missing data for network with index %d in %s\n",
                  logg.last_index, logf.file_name);
        }
        else
        { 
          s->sigma_block = logg.data['S'];
          w->param_block = logg.data['W'];
  
          net_setup_sigma_pointers (s, a);
          net_setup_param_pointers (w, a);

          if (op_N)
          { for (l = 0; l<a->N_layers; l++)
            { if (a->has_ho[l] && !keep[l])
              { for (k = 0; k<a->N_hidden[l]*a->N_outputs; k++)
                { w->ho[l][k] = 0;
                }
              }
            }
            if (a->has_io)
            { for (k = 0; k<a->N_inputs*a->N_outputs; k++)
              { w->io[k] = 0;
              }
            }
            if (a->has_bo)
            { for (k = 0; k<a->N_outputs; k++)
              { w->bo[k] = 0;
              }
            }
          }
    
          if (op_d && N_nets>=Max_median_nets)
          { fprintf(stderr,
              "Too many networks being used to find median (max %d)\n",
              Max_median_nets);
            exit(1);
          }

          for (i = 0; i<N_test; i++)
          { 
            net_func (&test_values[i], 0, a, w);
    
            if (have_targets) 
            { 
              net_model_prob (&test_values[i], 
                              test_targets + data_spec->N_targets * i, 
                              &curr_log_prob,
                              0, a, p, s, 0);
              if (op_r)
              { for (j = 0; j<data_spec->N_targets; j++)
                { tr = &data_spec->target_trans[j];
                  curr_log_prob += log(tr->scale);
                  if (tr->take_log)
                  { curr_log_prob -= log (data_inv_trans 
                      (test_targets[data_spec->N_targets*i+j], *tr));
                  }
                }
              }
    
              log_prob[i] = N_nets==0 ? curr_log_prob 
                                      : addlogs(log_prob[i],curr_log_prob); 
            }
    
            net_model_guess (&test_values[i], curr_targets, a, p, s, 0);
    
            for (j = 0; j<M_targets; j++)
            { 
              if (op_r) 
              { 
                tr = &data_spec->target_trans[j];

                curr_targets[j] = data_inv_trans(curr_targets[j], *tr);

                if (tr->take_log)
                {
                  if (op_n && a->data_model=='R' && (p->noise.alpha[0]!=0 
                       || p->noise.alpha[1]!=0 || p->noise.alpha[2]!=0))
                  { fprintf(stderr,"Predictive mean is undefined\n");
                    exit(1);
                  }

                  curr_targets[j] 
                    *= exp ((1/p->noise.width) / (2*tr->scale*tr->scale));
                }
              }

              sum_targets[M_targets*i+j] += curr_targets[j];
            }

            if (op_d)
            {
              rand_seed(101*logg.last_index+i);

              for (k = 0; k<Median_sample; k++)
              {
                net_model_guess (&test_values[i], curr_targets, a, p, s, 1);

                for (j = 0; j<M_targets; j++)
                { if (op_r)
                  { curr_targets[j] = data_inv_trans (curr_targets[j], 
                                                data_spec->target_trans[j]);
                  }
                  median_sample[i][j][ms_count+k] = curr_targets[j];
                }

              }
            }
          }

          ms_count += Median_sample;          
          N_nets += 1;
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
  
      logg.index['W'] = -1;  /* So they won't be mistaken for records from */
      logg.index['S'] = -1;  /* the new log file.                          */
    }
  
    if (N_nets==0) 
    { fprintf(stderr,"None of the networks specified were found\n");
      exit(1);
    }
  }

  /* Print headings. */

  if (!op_b) printf("\nNumber of networks used: %d\n",N_nets);

  if (!op_b && !op_a)
  { 
    printf("\nCase");

    if (op_i)
    { printf(" %-*s",7*a->N_inputs," Inputs");
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
      { printf(" %-*s",7*M_targets," Error^2");
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
 
    for (j = 0; j<data_spec->N_targets; j++) 
    { sq_error[j] = abs_error[j] = wrong[j] = 0;
      sq_error_sq[j] = abs_error_sq[j] = wrong_sq[j] = 0;
    }

    tsq_error = tabs_error =  0;
    tsq_error_sq = tabs_error_sq = 0;
  }

  for (i = 0; i<N_test; i++)
  { 
    if (put_in_blanks)
    { if (i>0 && test_values[i].i[0]!=test_values[i-1].i[0])
      { printf("\n");
      }
    }

    if (!op_a && !op_b) printf("%4d",i+1);
 
    if (!op_a && op_i)
    { printf(" ");
      for (j = 0; j<a->N_inputs; j++)
      { printf(" %6.2lf",test_values[i].i[j]);
      }
    }

    if (!op_a && op_t)
    { printf(" ");
      for (j = 0; j<data_spec->N_targets; j++)
      { double val;
        val = test_targets[data_spec->N_targets*i+j];
        if (op_r) val = data_inv_trans(val,data_spec->target_trans[j]);
        printf(" %6.2lf",val);
      }
      if (data_spec->N_targets==1) printf(" ");
    }

    if (op_p)
    { log_prob[i] -= log((double)N_nets);
      ave_log_prob += log_prob[i];
      ave_log_prob_sq += log_prob[i]*log_prob[i];
      if (!op_a) printf(" %9.3lf",log_prob[i]);
    }

    if (op_m)
    { 
      if (a->data_model=='C')
      { int g;
        g = 0;
        for (j = 1+op_z; j<a->N_outputs; j++)
        { if (sum_targets[M_targets*i+j] > sum_targets[M_targets*i+g]
           && (!op_z || g!=0 
                || sum_targets[M_targets*i+j] >
                   sum_targets[M_targets*i+0] + sum_targets[M_targets*i+1]))
          { g = j;
          }
        }
        guess[0] = g;
        if (have_targets)
        { error[0] = guess[0] != test_targets[data_spec->N_targets*i]
           && (!op_z || guess[0]!=0 || test_targets[data_spec->N_targets*i]!=1);
          wrong[0] += error[0];
          wrong_sq[0] += error[0]*error[0];
        }
      }

      if (a->data_model=='B')
      { for (j = 0; j<a->N_outputs; j++)
        { guess[j] = sum_targets[M_targets*i+j]/N_nets > 0.5;
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
        { printf(" %2.0lf",guess[j]);
        }
        if (data_spec->N_targets<3) printf("  ");
        if (have_targets)
        { printf(" ");
          if (data_spec->N_targets<2) printf("   ");
          for (j = 0; j<data_spec->N_targets; j++) 
          { printf(" %1.0lf",error[j]);
          }
        }
      }
    }

    if (op_n)
    { 
      double tot_error;

      tot_error = 0;

      for (j = 0; j<data_spec->N_targets; j++)
      { guess[j] = sum_targets[M_targets*i+j] / N_nets;
        if (have_targets)
        { double val, diff;
          val = test_targets[data_spec->N_targets*i+j];
          if (op_r) val = data_inv_trans(val,data_spec->target_trans[j]);
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
        for (j = 0; j<data_spec->N_targets; j++) 
        { printf(" %6.2lf",guess[j]);
        }
        if (have_targets)
        { printf(" ");
          for (j = 0; j<data_spec->N_targets; j++) 
          { printf(" %6.4lf",error[j]);
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
        { printf(" %6.2lf",guess[j]);
        }
        if (have_targets)
        { printf(" ");
          for (j = 0; j<data_spec->N_targets; j++) 
          { printf(" %6.4lf",error[j]);
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
      printf("Average log probability of targets: %9.3lf", a);
      if (N_test>1) 
      { printf("+-%.3lf", sqrt((ave_log_prob_sq/N_test-a*a) / (N_test-1)));
      }
      printf("\n\n");
    }

    if (op_m) 
    { double a;
      printf("Fraction of guesses that were wrong: ");
      for (j = 0; j<data_spec->N_targets; j++)
      { a = wrong[j] / N_test;
        printf(" %6.4lf",a);
        if (N_test>1) 
        { printf("+-%6.4lf", sqrt((wrong_sq[j]/N_test-a*a) / (N_test-1)));
        }
      }
      printf("\n");
    }

    if (op_n) 
    { 
      double a;

      printf("Average squared error guessing mean: ");

      for (j = 0; j<data_spec->N_targets; j++)
      { a = sq_error[j]/N_test;
        printf (" %8.5lf", a);
        if (N_test>1)
        { printf ("+-%.5lf", sqrt((sq_error_sq[j]/N_test - a*a) / (N_test-1)));
        }
      }
      printf("\n");

      if (data_spec->N_targets>1) 
      { a = tsq_error / N_test;
        printf("                                            (total %.5lf",a);
        if (N_test>1) 
        { printf("+-%.5lf",sqrt((tsq_error_sq/N_test - a*a) / (N_test-1)));
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
        printf(" %8.5lf", a);
        if (N_test>1)
        { printf("+-%.5lf", sqrt((abs_error_sq[j]/N_test - a*a) / (N_test-1)));
        }
      }
      printf("\n");

      if (data_spec->N_targets>1)
      { a = tabs_error / N_test;
        printf("                                            (total %.5lf",a);
        if (N_test>1) 
        { printf("+-%.5lf",sqrt((tabs_error_sq/N_test - a*a) / (N_test-1)));
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

static int dbl_cmp (const void *a, const void *b) 
{ return *(float*)a > *(float*)b ? +1 : *(float*)a < *(float*)b ? -1 : 0;
}

static float find_median
( float *data,		/* Array of data values */
  int n			/* Number of values in array */
)
{
  qsort (data, n, sizeof *data, dbl_cmp);
 
  return n%2==1 ? data[(int)(n/2)] 
                : (data[(int)(n/2)] + data[(int)((n-1)/2)]) / 2;
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf(stderr,
   "Usage: net-pred options { log-file range } [ / test-inputs [ test-targets ] ]\n");
  fprintf(stderr,
    "  Opt: i=inputs t=targets p=prob m=mode n=mean d=median\n");
  fprintf(stderr,
    "       r = raw b/B=bare a=aveonly 0-9=additive components\n");
  exit(1);
}
