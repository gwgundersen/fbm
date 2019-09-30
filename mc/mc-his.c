/* MC-HIS.C - Skeleton of program for Hamiltonian importance sampling. */

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


#define ECHO_ARGS 1	/* Set to 1 for debugging, normally 0 */
#define DIRECT_SUM 1	/* Set to 1 to disable fast summation, normally 0 */


/* CLOCKS_PER_SEC should be defined by <time.h>, but it seems that it
   isn't always.  1000000 seems to be the best guess for Unix systems. */

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000	/* Best guess */
#endif


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  static mc_iter it0, *it;
  static mc_traj tj0, *tj;
  static mc_temp_state temp_state;

  mc_dynamic_state ds;
  mc_value *start_q, *start_p, *save_q, *save_p;

  int n_traj, min_steps, max_steps;
  double stepsize, inv_temp, decay, mix;
  int modulus, pre_save, steps;

  log_file logf;
  log_gobbled logg;

  char **ap;
  char junk;

  double level, E0, E1, H0, H1;
  double *kinetic;

  int index;
  int i, j, k;

  unsigned old_clock; /* Theoretically, these should be of type clock_t, but  */
  unsigned new_clock; /* that type is inexplicably declared signed on most    */
                      /* systems, cutting the already-too-small range in half */

  /* Look at program arguments. */

  pre_save = 0;
  modulus = 1;
  steps = 1;
  mix = 0;

  if (argc<2) usage();

  logf.file_name = argv[1];
  
  ap = argv+2;

  if (*ap==0 || sscanf(*ap++,"%d%c",&n_traj,&junk)!=1 || n_traj<1
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

  if (*ap==0 || sscanf(*ap++,"%lf%c",&stepsize,&junk)!=1 || stepsize==0)
  { usage();
  }

  if (*ap!=0 && strcmp(*ap,"/")!=0)
  {
    if (sscanf(*ap++,"%d%c",&steps,&junk)!=1 || steps<1)
    { usage();
    }
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

  if (*ap!=0)
  { if (sscanf(*ap++,"%lf%c",&mix,&junk)!=1 || mix<0)
    { usage();
    }
  }

  if (*ap!=0) usage();

  if (ECHO_ARGS)  /* For debugging */
  { 
    printf("n-traj=%d min-steps=%d max-steps=%d modulus=%d pre_save=%d\n",
            n_traj, min_steps, max_steps, modulus, pre_save);
    printf("stepsize=%.3f steps=%d inv-temp=%.3f decay=%.3f mix=%.3f\n",
            stepsize, steps, inv_temp, decay, mix);
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

  tj = logg.data['t'];

  if (tj==0)
  { tj = &tj0;
    tj->type = 'L';
    tj->halfp = 1;
    tj->N_approx = 1;
  }

  it = logg.data['i'];

  if (it==0) it = &it0;

  if (it->approx_order[0]==0)
  { int na;
    na = tj->N_approx>0 ? tj->N_approx : -tj->N_approx;
    if (na>Max_approx) na = Max_approx;
    for (j = 0; j<na; j++) it->approx_order[j] = j+1;
  }

  ds.p = logg.data['p'];

  if (ds.p!=0)
  { if (logg.actual_size['p'] != ds.dim * sizeof (mc_value))
    { fprintf(stderr,"Momentum record has wrong size\n");
      exit(1);
    }
  }
  else
  { ds.p = chk_alloc (ds.dim, sizeof (mc_value));
  }

  /* Initialize for simulating dynamics. */

  ds.grad = chk_alloc (ds.dim, sizeof (mc_value));
  start_q = chk_alloc (ds.dim, sizeof (mc_value));
  start_p = chk_alloc (ds.dim, sizeof (mc_value));
  save_q  = chk_alloc (ds.dim, sizeof (mc_value));
  save_p  = chk_alloc (ds.dim, sizeof (mc_value));

  kinetic = chk_alloc (2*(max_steps-min_steps)+1, sizeof *kinetic);
  kinetic += max_steps-min_steps;

  if (stepsize>0)
  { mc_app_stepsizes(&ds);
    it->stepsize_factor = stepsize;
  }
  else
  { for (k = 0; k<ds.dim; k++) ds.stepsize[k] = 1;
    it->stepsize_factor = -stepsize;
  }

  mc_traj_init(tj,it);

  ds.temp_state = &temp_state;

  old_clock = clock();      /* start clock */
  it->time = 0;

  /* Simulate the requested number of trajectories. */

  for (i = 0; i<n_traj; i++)
  {
    /* Randomly generate the start state and save it. */

    mc_heatbath (&ds, 1/inv_temp, 0.0);
    if (!mc_app_zero_gen(&ds))
    { fprintf(stderr,
        "Application doesn't support Hamiltonian importance sampling\n");
      exit(1);
    }
    mc_value_copy (start_q, ds.q, ds.dim);
    mc_value_copy (start_p, ds.p, ds.dim);
    ds.know_grad    = 0;
    ds.know_pot     = 0;
    ds.know_kinetic = 0;

    /* Generate each of the states in turn, in the following order:
 
         -1, -2, ..., -(max_steps-min_steps), 0, 1, 2, ..., max_steps

       Note that if max_steps==min_steps, the sequence starts at 0. */

    it->move_point = 0;
    it->rejects = 0;
    it->proposals = 0;
    it->delta = 0;
    it->decay = decay;
    it->stepsize_factor = stepsize;

    j = max_steps>min_steps ? -1 : 0;

    while (j!=max_steps+1)
    {
      /* Create the next state. */

      if (j==0)
      { mc_value_copy (ds.q, start_q, ds.dim);
        for (k = 0; k<ds.dim; k++)
        { ds.p[k] = -start_p[k];
        }
        ds.know_grad    = 0;
        ds.know_pot     = 0;
        ds.know_kinetic = 0;
      }
      else 
      { 
        if (j<0)
        { for (k = 0; k<ds.dim; k++)
          { ds.p[k] /= decay;
          }
          ds.know_kinetic = 0;
          if (mix!=0)
          { mc_mix_momentum (&ds, mix); 
          }
        }

        ds.temp_state->inv_temp = 0;
        mc_app_energy(&ds,1,1,&E0,0);
        level = E0 + rand_exp();
        mc_value_copy (save_q, ds.q, ds.dim);
        mc_value_copy (save_p, ds.p, ds.dim);

        ds.temp_state->inv_temp = -1;
        mc_app_energy(&ds,1,1,&H0,0);
        H0 += mc_kinetic_energy(&ds);

        mc_trajectory (&ds, steps, 0);

        mc_app_energy(&ds,1,1,&H1,0);
        H1 += mc_kinetic_energy(&ds);

        ds.temp_state->inv_temp = 0;
        mc_app_energy(&ds,1,1,&E1,0);

        if (E1>level || (H1-H0)>1000 || (H0-H1)>1000)
        { mc_value_copy (ds.q, save_q, ds.dim);
          for (k = 0; k<ds.dim; k++)
          { ds.p[k] = -save_p[k];
          }
          it->rejects += 1;
        }
        it->proposals += 1;
        it->delta = E1-E0;

        if (j>0)
        { if (mix!=0)
          { mc_mix_momentum (&ds, mix); 
          }
          for (k = 0; k<ds.dim; k++)
          { ds.p[k] *= decay;
          }
          ds.know_kinetic = 0;
        }
      }

      /* Save kinetic energy for possible start states. */

      if (j<=max_steps-min_steps)
      { ds.kinetic_energy = mc_kinetic_energy(&ds);
        kinetic[j] = ds.kinetic_energy;
        ds.know_kinetic = 1;
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
          it->log_weight = - inv_temp * kinetic[j-min_steps];
          for (k = 1; k<=max_steps-min_steps; k++)
          { it->log_weight = addlogs (it->log_weight,
              - inv_temp * kinetic[j-min_steps-k] - k*ds.dim*log(decay));
          } 
          it->log_weight = - it->log_weight
                           - 0.5*ds.dim*log(inv_temp)
                           + min_steps*ds.dim*log(decay)
                           + log(max_steps-min_steps+1.0);
  
          ds.temp_state->inv_temp = -1;
          mc_app_energy(&ds,1,1,&ds.pot_energy,0);
          ds.kinetic_energy = mc_kinetic_energy(&ds);
          it->log_weight -= ds.pot_energy + ds.kinetic_energy;
        }
  
        /* Record position along trajectory for this state. */
  
        it->move_point = j;

        /* See how long it took to create this state. */
    
        new_clock = clock(); 
        it->time += (int) (0.5 + 
               (1000.0 * (unsigned) (new_clock - old_clock)) / CLOCKS_PER_SEC);
        old_clock = new_clock;

        /* Write records to log file. */

        mc_app_save(&ds,&logf,index);
  
        logf.header.type = 'p';
        logf.header.index = index;
        logf.header.size = ds.dim * sizeof (mc_value);
        log_file_append (&logf, ds.p);
  
        logf.header.type = 'r';
        logf.header.index = index;
        logf.header.size = sizeof (rand_state);
        log_file_append (&logf, rand_get_state());
  
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
   "Usage: xxx-his log-file n-traj min-steps max-steps [ [-]modulus ]\n");
  fprintf (stderr,
   "                      / stepsize [ steps ] / inv-temp decay [ mix ]\n");
  exit(1);
}
