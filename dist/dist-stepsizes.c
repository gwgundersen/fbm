/* DIST-STEPSIZES.C - Specialized stepsize program for 'dist' module. */

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
#include "log.h"
#include "mc.h"
#include "formula.h"
#include "dist.h"

static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  mc_dynamic_state ds;
  dist_spec *dst;

  log_file logf;
  log_gobbled logg;

  mc_value *grad, *bgrad, *fgrad;
  mc_value *seconds, *stepsizes;

  int show, set, index, any_specs;
  double delta; 
 
  int a, c, i, n;
  char *e, *p;

  /* Look at program arguments. */

  show = 0;
  set = 0;

  if (argc<2) usage();

  logf.file_name = argv[1];

  if (argc>2 && strchr(argv[2],'=')!=0)
  { set = 1;
  }
  else
  { index = -1;
    if (argc>2)
    { if ((index = atoi(argv[2]))<0 || index==0 && strcmp(argv[2],"0")!=0) 
      { usage();
      }
    }
    if (argc<=3)
    { show = 1;
    }
    else
    { if (argc!=4 || (delta = atof(argv[3]))<=0) usage();
    }
  }

  /* Open log file and read records with indexes less than zero. */

  log_file_open (&logf, set);

  log_gobble_init(&logg,0);
  mc_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  /* Do generic display of actual and optimal stepsizes if requested. */

  if (!set && !show)
  { 
    /* Read records with desired index. */
  
    while (!logf.at_end && logf.header.index!=index)
    { log_gobble(&logf,&logg);
    }
  
    if (logf.at_end) 
    { fprintf(stderr,"No records with that index in log file\n");
      exit(1);
    }
  
    log_gobble(&logf,&logg);
  
    /* Initialize using records found. */
  
    mc_app_initialize(&logg,&ds);

    /* Allocate space. */

    grad = chk_alloc (ds.dim, sizeof *grad);
    bgrad = chk_alloc (ds.dim, sizeof *bgrad);
    fgrad = chk_alloc (ds.dim, sizeof *fgrad);
    seconds = chk_alloc (ds.dim, sizeof *seconds);

    /* Get stepsizes selected by application. */

    mc_app_stepsizes(&ds);

    /* Compute second derivatives to get comparison stepsizes. */

    for (i = 0; i<ds.dim; i++)
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

    /* Print results. */

    printf(
     "\nActual stepsizes and stepsizes based on second derivatives\n\n");
    printf("Log file: %s  Index: %d  Delta: %.6f\n\n",
            logf.file_name, index, delta);
    printf("  Coord  Stepsize   Back. Grad.  Forw. Grad.   Stepsize\n\n");

    for (i = 0; i<ds.dim; i++)
    { 
      if (seconds[i]>0)
      { printf ("%6d %+10.5f   %+11.4e  %+11.4e   %+10.5f\n",
                i, ds.stepsize[i], bgrad[i], fgrad[i], 1/sqrt(seconds[i]));
      }
      else
      { printf ("%6d %+10.5f   %+11.4e  %+11.4e   **********\n",
                i, ds.stepsize[i], bgrad[i], fgrad[i]);
      }
    }
  
    printf("\n");
  }

  /* Look at specification if we're going to set or show stepsizes by name. */

  if (show || set)
  { 
    dst = logg.data['d'];
  
    if (dst==0)
    { fprintf(stderr,"No distribution specification in log file\n");
      exit(1);
    }

    (void) formula (dst->energy, 1, 0, 0);
    if (dst->Bayesian)
    { (void) formula (dst->energy + strlen(dst->energy) + 1, 1, 0, 0);
    }
  }

  /* Show stepsizes specified if requested. */

  if (show)
  {
    any_specs = 0;

    if (index==-1)
    {
      while (!logf.at_end)
      { log_gobble(&logf,&logg);
      }
  
      if (logg.data['s']!=0)
      { printf("\nLAST STEPSIZE SPECIFICATIONS IN LOG FILE:\n\n");
        any_specs = 1;
      }
  
    }
    else
    {
      while (!logf.at_end && logf.header.index<=index)
      { log_gobble(&logf,&logg);
      }

      if (logg.data['s']!=0 && logg.index['s']==index)
      { printf("\nSTEPSIZE SPECIFICATIONS AT INDEX %d IN LOG FILE:\n\n",index);
        any_specs = 1;
      }
    }

    if (!any_specs)
    { printf("\nNo stepsize specification found\n\n");
    }
    else
    { 
      if (logg.actual_size['s']!=dist_count_vars()*sizeof(double))
      { fprintf(stderr,"Stepsize record is the wrong size!\n");
        exit(1);
      }

      dist_unpack_vars(logg.data['s']);
      dist_print_vars();
      printf("\n");
    }
  }

  /* Set stepsizes to specified values if requested. The values are not 
     allowed to be general formulas, since it's hard to arrange for 
     variable references to be prohibited. */

  if (set)
  {
    n = dist_count_vars();
    stepsizes = chk_alloc (n, sizeof(double));

    for (i = 0; i<n; i++) 
    { stepsizes[i] = 1.0;
    }

    dist_unpack_vars(stepsizes);

    for (a = 2; argv[a]; a++)
    { 
      if (!strchr(argv[a],'='))
      { usage();
      }
  
      e = formula_def(argv[a],&c,&i);
      
      if (!formula_var_exists[c][i])
      { *strchr(argv[a],'=') = 0;
        fprintf(stderr,"Not a specified state variable: %s\n",argv[a]);
        exit(1);
      }
  
      formula_var[c][i] = strtod(e,&p);
      if (*p!=0) usage();

      if (formula_var[c][i]<=0) 
      { fprintf(stderr,"Stepsizes must be positive\n");
        exit(1);
      }
    }

    dist_pack_vars(stepsizes);

    logf.header.type = 's';
    logf.header.index = logf.at_end || logf.header.index<0 ? 0 
                         : logf.header.index;
    logf.header.size = n * sizeof(double);
    log_file_append(&logf,stepsizes);
  }

  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

static void usage(void)
{ fprintf(stderr,
    "Usage: dist-stepsizes log-file index delta\n");
  fprintf(stderr,
    "   or: dist-stepsizes log-file { var=value }\n");
  fprintf(stderr,
    "   or: dist-stepsizes log-file [ index ] (to display stored stepsizes)\n");
  exit(1);
}
