/* MC-MIS.C - Skeleton of program for Metropolis importance sampling. */

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
#include <time.h>

#include "misc.h"
#include "rand.h"
#include "log.h"
#include "mc.h"


#define ECHO_ARGS 1		/* Set to 1 for debugging, normally 0 */

#define Max_stepsizes 10	/* Maximum number of stepsize/count args */

/* CLOCKS_PER_SEC should be defined by <time.h>, but it seems that it
   isn't always.  1000000 seems to be the best guess for Unix systems. */

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000	/* Best guess */
#endif


static void usage(void);


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  static mc_iter it0, *it;
  static mc_therm_state th0, *th;
  static mc_temp_state temp_state;

  mc_dynamic_state ds;
  double *start_q, *save_q;
  double start_thm;

  int n_runs, min_steps, max_steps, blocksize;
  double stepsize[Max_stepsizes];
  int count[Max_stepsizes];
  int n_stepsizes;
  double inv_temp, decay;
  int modulus, pre_save;

  log_file logf;
  log_gobbled logg;

  char **ap;
  char junk;

  double ipr, E0, E1, H0, H1;
  double *thm;

  int index;
  int i, j, k, s, c, b;

  unsigned old_clock; /* Theoretically, these should be of type clock_t, but  */
  unsigned new_clock; /* that type is inexplicably declared signed on most    */
                      /* systems, cutting the already-too-small range in half */

  /* Look at program arguments. */

  pre_save = 0;
  modulus = 1;

  if (argc<2) usage();

  logf.file_name = argv[1];
  
  ap = argv+2;

  if (*ap==0 || sscanf(*ap++,"%d%c",&n_runs,&junk)!=1 || n_runs<1
   || *ap==0 || sscanf(*ap++,"%d%c",&min_steps,&junk)!=1 || min_steps<0
   || *ap==0 || sscanf(*ap++,"%d%c",&max_steps,&junk)!=1 || max_steps<min_steps)
  { usage();
  }

  if (*ap!=0 && strchr("-0123456789",**ap)!=0)
  { if (sscanf(*ap++,"%d%c",&modulus,&junk)!=1 || modulus==0)
    { usage();
    }
  }

  if (modulus<0)
  { pre_save = 1;
    modulus = -modulus;
  }

  if (*ap==0 || strcmp(*ap++,"/")!=0)
  { usage();
  }

  if (*ap!=0 && strcmp(*ap,"all")==0)
  { blocksize = 1000000000;
    ap += 1;
  }
  else if (*ap==0 || sscanf(*ap++,"%d%c",&blocksize,&junk)!=1 || blocksize<=0)
  { usage();
  }

  n_stepsizes = 0;
  while (*ap!=0 && strcmp(*ap,"/")!=0) 
  {
    if (n_stepsizes>=Max_stepsizes)
    { fprintf(stderr,"Too many stepsize/count specifications (max %d)\n",
        Max_stepsizes);
      exit(1);
    }

    if (sscanf(*ap++,"%lf%c",&stepsize[n_stepsizes],&junk)!=1 
     || stepsize[n_stepsizes]<=0)
    { usage();
    }
    
    if (*ap==0 || strcmp(*ap,"/")==0) 
    { count[n_stepsizes] = 1;
    }
    else if (sscanf(*ap++,"%d%c",&count[n_stepsizes],&junk)!=1 
          || count[n_stepsizes]<=0)
    { usage();
    }

    n_stepsizes += 1;
  }
 
  if (n_stepsizes==0)
  { usage();
  }

  if (*ap==0 || strcmp(*ap++,"/")!=0)
  { usage();
  }

  if (*ap==0 || **ap=='/' && (sscanf(1+*ap++,"%lf%c",&inv_temp,&junk)!=1
                               || inv_temp<=0 || (inv_temp=1/inv_temp)<1)
             || **ap!='/' && (sscanf(*ap++,"%lf%c",&inv_temp,&junk)!=1
                               || inv_temp<=0 || inv_temp>1)
   || *ap==0 || sscanf(*ap++,"%lf%c",&decay,&junk)!=1 || decay<=0 || decay>1)
  { usage();
  }

  if (*ap!=0) usage();

  if (ECHO_ARGS)  /* For debugging */
  { 
    printf("n-runs=%d min-steps=%d max-steps=%d modulus=%d pre_save=%d\n",
            n_runs, min_steps, max_steps, modulus, pre_save);
    printf("blocksize=%d  stepsize/count:",blocksize);
    for (i = 0; i<n_stepsizes; i++)
    { printf("  %.3f %d",stepsize[i],count[i]);
    }
    printf("\n");
    printf("inv-temp=%.3f decay=%.3f\n",inv_temp, decay);
    fflush(stdout);
  }

  /* Open log file and read all records. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  mc_record_sizes(&logg);

  while (!logf.at_end)
  { log_gobble(&logf,&logg);
  }

  index = log_gobble_last(&logf,&logg);

  /* Look at what records we have.  Use them if we have them; use defaults
     if we don't. */

  ds.aux_dim = 0;
  ds.aux = 0;

  mc_app_initialize(&logg,&ds);

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  it = logg.data['i'];
  if (it==0) it = &it0;

  th = logg.data['h'];
  if (th==0) th = &th0;

  /* Initialize for simulation. */

  start_q = chk_alloc (ds.dim, sizeof (mc_value));
  save_q  = chk_alloc (ds.dim, sizeof (mc_value));

  thm = chk_alloc (2*(max_steps-min_steps)+1, sizeof *thm);
  thm += max_steps-min_steps;

  mc_app_stepsizes(&ds);

  ds.temp_state = &temp_state;

  old_clock = clock();      /* start clock */
  it->time = 0;

  /* Simulate the requested number of runs. */

  for (i = 0; i<n_runs; i++)
  {
    /* Randomly generate the start state and save it. */

    if (!mc_app_zero_gen(&ds))
    { fprintf(stderr,
        "Application doesn't support Metropolis importance sampling\n");
      exit(1);
    }
    mc_value_copy (start_q, ds.q, ds.dim);
    ds.know_pot = 0;

    th->tq = rand_exp() / inv_temp;
    start_thm = th->tq;

    /* Generate each of the states in turn, in the following order:
 
         -1, -2, ..., -(max_steps-min_steps), 0, 1, 2, ..., max_steps

       Note that if max_steps==min_steps, the sequence starts at 0. */

    it->move_point = 0;
    it->rejects = 0;
    it->proposals = 0;
    it->delta = 0;
    it->decay = decay;

    j = max_steps>min_steps ? -1 : 0;

    while (j!=max_steps+1)
    {
      /* Create the next state. */

      if (j==0)
      { mc_value_copy (ds.q, start_q, ds.dim);
        ds.know_pot = 0;
        th->tq = start_thm;
      }
      else 
      { 
        if (j<0)
        { th->tq /= decay;
        }

        for (s = 0; s<n_stepsizes; s++)
        {
          for (c = 0; c<count[s]; c++)
          { 
            for (b = 0; b<ds.dim; b+=blocksize)
            {
              mc_value_copy (save_q, ds.q, ds.dim);

              ds.temp_state->inv_temp = 0;
              mc_app_energy(&ds,1,1,&E0,0);
              ds.temp_state->inv_temp = -1;
              mc_app_energy(&ds,1,1,&H0,0);

              for (k = 0; k<blocksize && b+k<ds.dim; k++)
              { ds.q[b+k] += rand_gaussian() * ds.stepsize[b+k] * stepsize[s];
              }

              mc_app_energy(&ds,1,1,&H1,0);
              ds.temp_state->inv_temp = 0;
              mc_app_energy(&ds,1,1,&E1,0);

              if (rand_uniform()>=exp(E0-E1) || H1-H0>=th->tq)
              { mc_value_copy (ds.q, save_q, ds.dim);
                it->rejects += 1;
              }
              else
              { th->tq -= H1-H0;
              }

              it->proposals += 1;
              it->delta = H1-H0;
            }
          }
        }

        if (j>0)
        { th->tq *= decay;
        }
      }

      /* Save thermostat value for possible start state. */

      if (j<=max_steps-min_steps)
      { thm[j] = th->tq;
      }

      /* Find weight for state and write it to the log file if appropriate. */

      if (j%modulus==0 && (pre_save || j>=min_steps))
      { 
        /* Compute weight for this state. */
  
        if (j<min_steps)
        { it->log_weight = -1e30; /* effectively zero */
        }
        else 
        { 
          for (k = min_steps; k<=max_steps; k++)
          { ipr = k==min_steps ? -inv_temp*thm[j-k] - k*log(decay)
                  : addlogs(ipr, -inv_temp*thm[j-k] - k*log(decay));
          } 
          ipr += log(inv_temp) - log(max_steps-min_steps+1.0);

          ds.temp_state->inv_temp = -1;
          mc_app_energy(&ds,1,1,&ds.pot_energy,0);
          it->log_weight = - (ds.pot_energy + th->tq) - ipr;
        }

        /* Record position within run for this state. */
  
        it->move_point = j;

        /* See how long it took to create this state. */
    
        new_clock = clock(); 
        it->time += (int) (0.5 + 
               (1000.0 * (unsigned) (new_clock - old_clock)) / CLOCKS_PER_SEC);
        old_clock = new_clock;

        /* Write records to log file. */

        mc_app_save(&ds,&logf,index);
  
        logf.header.type = 'r';
        logf.header.index = index;
        logf.header.size = sizeof (rand_state);
        log_file_append (&logf, rand_get_state());
  
        logf.header.type = 'h';
        logf.header.index = index;
        logf.header.size = sizeof *th;
        log_file_append (&logf, th);
  
        logf.header.type = 'i';
        logf.header.index = index;
        logf.header.size = sizeof *it;
        log_file_append (&logf, it);

        it->rejects = 0;
        it->proposals = 0;

        index += 1;
      }

      /* Go on to the next state. */

      if (j>=0)
      { j += 1;
      }
      else if (j==-(max_steps-min_steps))
      { j = 0;
      }
      else
      { j -= 1;
      }
    }
  }

  log_file_close(&logf);

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf (stderr, 
 "Usage: xxx-mis log-file n-runs min-steps max-steps [ [-]modulus ]\n");
  fprintf (stderr,
 "         / blocksize { stepsize count } stepsize [ count ] / inv-temp decay\n"
  );
  exit(1);
}
