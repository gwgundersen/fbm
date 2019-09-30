/* MC-STEPSIZES.C - Skeleton program to display and evaluate stepsizes. */

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


/* MAIN PROGRAM. */

void main
( int argc,
  char **argv
)
{
  mc_dynamic_state ds;

  log_file logf;
  log_gobbled logg;

  mc_value *grad, *bgrad, *fgrad;
  mc_value *seconds;

  double delta; 
  int index;
  int dosec;
 
  int i;

  /* Look at program arguments. */

  dosec = argc==4;

  if (argc!=3 && argc!=4
   || (index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0
   || dosec && (delta = atof(argv[3]))<=0)
  { fprintf(stderr,"Usage: xxx-stepsizes log-file index [ delta ]\n");
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
  bgrad = chk_alloc (ds.dim, sizeof *bgrad);
  fgrad = chk_alloc (ds.dim, sizeof *fgrad);

  if (dosec) seconds = chk_alloc (ds.dim, sizeof *seconds);

  /* Get stepsizes selected by application. */

  mc_app_stepsizes(&ds);

  /* Compute second derivatives to get comparison stepsizes. */

  if (dosec)
  { for (i = 0; i<ds.dim; i++)
    { double t;
      t = ds.q[i];
      ds.q[i] = t - delta/2;
      mc_app_energy (&ds, 1, 1, (double *) 0, grad);
      bgrad[i] = grad[i];
      ds.q[i] = t + delta/2;
      mc_app_energy (&ds, 1, 1, (double *) 0, grad);
      fgrad[i] = grad[i];
      ds.q[i] = t;
      seconds[i] = (fgrad[i] - bgrad[i]) / delta;
    }
  }

  /* Print results. */

  if (dosec)
  { printf(
     "\nStepsizes selected by application and from second derivatives\n\n");
    printf("Log file: %s  Index: %d  Delta: %.6lf\n\n",
            logf.file_name, index, delta);
    printf("  Coord  Stepsize   Back. Grad.  Forw. Grad.   Stepsize\n\n");
  }
  else
  { printf("\nStepsizes selected by application\n\n");
    printf("Log file: %s  Index: %d\n\n", logf.file_name, index);
    printf("  Coord  Stepsize\n\n");
  }

  for (i = 0; i<ds.dim; i++)
  { 
    if (dosec)
    { if (seconds[i]>0)
      { printf ("%6d %+10.5lf   %+11.4le  %+11.4le   %+10.5lf\n",
                i, ds.stepsize[i], bgrad[i], fgrad[i], 1/sqrt(seconds[i]));
      }
      else
      { printf ("%6d %+10.5lf   %+11.4le  %+11.4le   **********\n",
                i, ds.stepsize[i], bgrad[i], fgrad[i]);
      }
    }
    else
    { printf ("%6d %+10.5lf\n", i, ds.stepsize[i]);
    }
  }

  printf("\n");

  exit(0);
}
