/* NET-APPROX.C - Specify quadratic approximation to the log likelihood. */

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
#include "prior.h"
#include "data.h"
#include "model.h"
#include "net.h"

static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  net_arch *a;
  net_flags *flgs;

  log_file logf;
  log_gobbled logg;

  char *approx_file;
  FILE *af;

  int value_count;
  double *approx;
  double junk;

  int i, j, n;

  /* Look at arguments. */
  
  if (argc!=3) usage();

  logf.file_name = argv[1];
  approx_file = argv[2];

  /* Open log file and read network architecture. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  net_record_sizes(&logg);

  if (!logf.at_end && logf.header.index==-1)
  { log_gobble(&logf,&logg);
  }

  a = logg.data['A'];
  flgs = logg.data['F'];

  if (a==0)
  { fprintf(stderr,"No architecture specification in log file\n");
    exit(1);
  }

  /* Allocate space for parameters of approximation. */

  n = net_setup_param_count(a,flgs);

  approx = chk_alloc (1+n+n*n, sizeof *approx);

  /* Open approximation file and read parameters of approximation from it. */

  af = fopen(approx_file,"r");
  if (af==NULL)
  { fprintf(stderr,
           "Can't open approximation file for reading (%s)\n",approx_file);
    exit(1);
  }

  for (i = 0; i<1+n+n*n; i++)
  { if (fscanf(af,"%lf",&approx[i])!=1)
    { fprintf(stderr,"Error reading from approximation file\n");
      exit(1);
    }
  }

  if (fscanf(af,"%lf",junk)!=EOF)
  { fprintf(stderr,"Excess data in approximation file\n");
    exit(1);
  }

  /* Write approximation record to log file. */

  logf.header.type = 'Q';
  logf.header.index = -1;
  logf.header.size = (1+n+n*n) * sizeof *approx;
  log_file_append (&logf, approx);

  log_file_close (&logf);
  
  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf (stderr, "Usage: net-approx log-file approximation-file\n");
  exit(1);
}
