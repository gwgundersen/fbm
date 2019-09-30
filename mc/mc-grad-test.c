/* MC-GRAD-TEST.C - Skeleton program to test energy gradient computation. */

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


/* THE FOLLOWING IS NEEED BY THE PLOTTING ROUTINES. */

enum { PLT, TBL, HIST } program_type;


/* MAIN PROGRAM. */

void main
( int argc,
  char **argv
)
{
  mc_dynamic_state ds;

  log_file logf;
  log_gobbled logg;

  mc_value *grad;

  double E, *EF, *EB;

  double delta; 
  int index;

  double maxlr;
  int negr;
 
  int i;

  /* Look at program arguments. */

  if (argc!=4 
   || (index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0
   || (delta = atof(argv[3]))<=0)
  { fprintf(stderr,"Usage: xxx-grad-test log-file index delta\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read records with indexes less than zero and with
     index equal to that at the end. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  mc_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  /* Read records with desired index. */

  while (!logf.at_end && logf.header.index!=index)
  { log_file_forward(&logf);
  }

  if (logf.at_end) 
  { fprintf(stderr,"No records with that index in log file\n");
    exit(1);
  }

  log_gobble(&logf,&logg);

  /* Initialize using records found. */

  mc_app_initialize(&logg,&ds);

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  grad = chk_alloc (ds.dim, sizeof *grad);

  EF = chk_alloc (ds.dim, sizeof *EF);
  EB = chk_alloc (ds.dim, sizeof *EB);

  /* Find gradient using application-specific gradient procedure. */

  mc_app_energy (&ds, 1, 1, &E, grad);

  /* Gather differences. */

  for (i = 0; i<ds.dim; i++)
  { double t;
    t = ds.q[i];
    ds.q[i] = t + delta/2;
    mc_app_energy (&ds, 1, 1, &EF[i], 0);
    ds.q[i] = t - delta/2;
    mc_app_energy (&ds, 1, 1, &EB[i], 0);
    ds.q[i] = t;
  }

  /* Print results. */

  printf("\nComparison of computed gradient with finite differences\n\n");

  printf("Log file: %s  Index: %d  Delta: %.6f  Energy: %.2f\n\n",
          logf.file_name, index, delta, E);

  printf(
   "  Coord  Value   Forw. Diff.  Back. Diff.   Diff. Grad.  Comp. Grad.   Log Ratio\n\n");

  negr = 0;
  maxlr = 0;

  for (i = 0; i<ds.dim; i++)
  { 
    double r, lr;

    printf ("%6d %+7.2f   %+11.4e  %+11.4e   %+11.4e  %+11.4e", 
            i, ds.q[i], EF[i]-E, E-EB[i], (EF[i]-EB[i])/delta, grad[i]);

    r = ((EF[i]-EB[i])/delta) / grad[i];
    if (r<=0) 
    { printf("   ********\n");
      negr = 1;
    }
    else
    { lr = log(r);
      printf("   %+8.4f\n", lr);
      if (lr<0) lr = -lr;
      if (lr>maxlr) maxlr = lr;
    }
  }

  if (negr) printf("\nLargest absolute log ratio: ********\n\n");
  else      printf("\nLargest absolute log ratio: %.4f\n\n",maxlr);

  exit(0);
}
