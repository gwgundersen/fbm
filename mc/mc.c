/* MC.C - Skeleton of program to run Markov chain Monte Carlo simulation. */

/* Copyright (c) 1995, 1996, 1998 by Radford M. Neal 
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
 * The features allowing the number of iterations to be specified in terms of
 * cpu time are adapted from modifications done by Carl Edward Rasmussen, 1995.
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
#include "quantities.h"


/* CLOCKS_PER_SEC should be defined by <time.h>, but it seems that it
   isn't always.  1000000 seems to be the best guess for Unix systems. */

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000	/* Best guess */
#endif


/* THE FOLLOWING IS NEEED BY THE PLOTTING ROUTINES. */

enum { PLT, TBL, HIST } program_type = PLT;


static void usage(void);


/* MAIN PROGRAM. */

void main
( int argc,
  char **argv
)
{
  int index, max, modulus;

  static mc_ops *ops;
  static mc_traj tj0, *tj;
  static mc_iter it0, *it;

  mc_temp_sched *sch;
  mc_dynamic_state ds;

  char *quantities;
  int N_quantities;

  quantities_described qd;

  double temperature, decay;

  log_file logf;
  log_gobbled logg;

  char **ap;

  int j, na;
  int timelimit;

  unsigned old_clock; /* Theoretically, these should be of type clock_t, but  */
  unsigned new_clock; /* that type is inexplicably declared signed on most    */
                      /* systems, cutting the already-too-small range in half */

  /* Look at program arguments. */

  timelimit = 0;
  N_quantities = 0;
  quantities = 0;
  modulus = 1;
  temperature = 1;
  decay = -1;

  if (argc<3) usage();

  ap = argv+2;

  if (**ap == '@') timelimit = 1; 

  if ((max = atoi(**ap=='@' ? *ap+1 : *ap)) <= 0) usage();
  ap += 1;

  if (*ap!=0 && strchr("0123456789",**ap)!=0 
       && ((modulus = atoi(*ap++))<=0)) usage();

  if (*ap!=0 && strcmp(*ap,"/")!=0) 
  { quantities = *ap++; 
  }

  if (*ap!=0 && strcmp(*ap,"/")==0)
  { ap += 1;
    if (*ap==0 || (decay = atof(*ap))<=0 && strcmp(*ap,"0")!=0) usage();
    ap += 1;
    if (*ap!=0)
    { if ((temperature = atof(*ap))<=0 && strcmp(*ap,"0")!=0) usage();
      ap += 1;
    }
  }

  if (*ap!=0) usage();

  logf.file_name = argv[1];

  /* Open log file and read all records. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  mc_record_sizes(&logg);

  while (!logf.at_end)
  { log_gobble(&logf,&logg);
  }

  index = log_gobble_last(&logf,&logg);

  if (!timelimit && index>max)
  { fprintf(stderr,"Iterations up to %d already exist in log file\n",max);
    exit(1);
  }

  /* Look at what records we have.  Use them if we have them; use defaults
     if we don't. */

  ds.aux_dim = 0;
  ds.aux = 0;

  mc_app_initialize(&logg,&ds);

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  ops = logg.data['o'];

  if (ops==0)
  { fprintf(stderr,"No Markov chain operations specified\n");
    exit(1);
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
  { na = tj->N_approx>0 ? tj->N_approx : -tj->N_approx;
    if (na>Max_approx) na = Max_approx;
    for (j = 0; j<na; j++) it->approx_order[j] = j+1;
  }

  sch = logg.data['m'];

  ds.p = logg.data['p'];

  if (ds.p!=0)
  { if (logg.actual_size['p'] != ds.dim * sizeof (mc_value))
    { fprintf(stderr,"Momentum record has wrong size (in mc.c)\n");
      exit(1);
    }
  }

  logg.req_size['p'] = ds.dim * sizeof (mc_value);

  ds.temp_state = logg.data['b'];

  if (ds.temp_state!=0)
  { ds.temp_index = mc_temp_index (sch, ds.temp_state->inv_temp);
  }

  /* Initialize for performing iterations. */

  ds.grad = 0;

  ds.know_grad    = 0;
  ds.know_pot     = 0;
  ds.know_kinetic = 0;

  mc_iter_init(&ds,ops,tj,sch);

  if (quantities)
  { 
    if (quantities[0]=='t' && quantities[1]=='t' && quantities[2]!=0 
     && strchr("12",quantities[2])!=0 && quantities[3]==0)
    { N_quantities = -(quantities[2]-'0');
    }
    else
    {
      static char *zero = 0, **null_argv = &zero;
      int a;

      for (a = 0; quant_app_arguments[a]; a++)
      { (*quant_app_arguments[a])(&null_argv);
      }

      for (a = 0; quant_app_initialize[a]; a++)
      { (*quant_app_initialize[a])(&logg);
      }

      N_quantities = quantities_requested(qd,quantities,1);
      quantities_available(qd,&logg);  
    }
  }

  it->temperature = temperature;
  it->decay = decay;

  it->proposals = 0;
  it->rejects = 0;
  it->slice_calls = 0;
  it->slice_evals = 0;

  old_clock = clock();      /* start clock */

  /* Perform Markov chain iterations. */

  for (; (!timelimit && index<=max) || (timelimit && it->time<60000*max); 
         index++)
  {
    mc_iteration (&ds, it, &logg, qd, N_quantities);

    new_clock = clock(); 
    it->time += (int) (0.5 + 
             (1000.0 * (unsigned) (new_clock - old_clock)) / CLOCKS_PER_SEC);
    old_clock = new_clock;

    if (index%modulus==0)
    { 
      mc_app_save(&ds,&logf,index);

      if (ds.p!=0)
      { logf.header.type = 'p';
        logf.header.index = index;
        logf.header.size = ds.dim * sizeof (mc_value);
        log_file_append (&logf, ds.p);
      }

      if (ds.temp_state!=0)
      { logf.header.type = 'b';
        logf.header.index = index;
        logf.header.size = sizeof (mc_temp_state);
        log_file_append (&logf, ds.temp_state);
      }

      logf.header.type = 'i';
      logf.header.index = index;
      logf.header.size = sizeof *it;
      log_file_append (&logf, it);

      logf.header.type = 'r';
      logf.header.index = index;
      logf.header.size = sizeof (rand_state);
      log_file_append (&logf, rand_get_state());

      it->proposals = 0;
      it->rejects = 0;
      it->slice_calls = 0; 
      it->slice_evals = 0;
    }
  }

  log_file_close(&logf);

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf (stderr, 
    "Usage: xxx-mc log-file [\"@\"]iteration [ save-mod ]\n"
    "              [ quantities | \"tt[12]\" ] [ / decay [ temperature ] ]\n");
  exit(1);
}
