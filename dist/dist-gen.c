/* DIST-GEN.C - Program to sample from prior distribution. */

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
#include "dist.h"
#include "formula.h"
#include "rand.h"


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  dist_spec *dst;
  char *a, *f;
  double *q;

  int max_index, index;
  int c, i, n;

  log_file logf;
  log_gobbled logg;
 
  /* Look at program arguments. */

  max_index = 0;

  if (argc<2 || argc>3
   || argc>2 && (max_index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0)
  { fprintf(stderr,"Usage: dist-gen log-file [ max-index ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  logg.req_size['d'] = sizeof *dst;
  logg.req_size['r'] = sizeof (rand_state);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  /* Look at specification, to find out what state variables exist. */

  if ((dst = logg.data['d'])==0)
  { fprintf(stderr,"No distribution specification in log file\n");
    exit(1);
  }

  if (!dst->Bayesian)
  { fprintf(stderr,
      "Generation from the prior can be done only with Bayesian models\n");
    exit(1);
  }

  a = dst->energy + strlen(dst->energy) + 1;
  if (dst->Bayesian) a += strlen(a) + 1;

  while (*a)
  { f = formula_def(a,&c,&i);
    formula_var[c][i] = formula(f,0,1,0);
    formula_var_exists[c][i] = 1;
    a += strlen(a) + 1;
  }

  (void) formula (dst->energy, 1, 0, 0);
  (void) formula (dst->energy + strlen(dst->energy) + 1, 1, 0, 0);

  /* Read last records in log file to see where to start, and to get random
     number state left after last network was generated. */

  index = log_gobble_last(&logf,&logg);

  if (logg.last_index<0) 
  { index = 0;
  }

  if (index>max_index)
  { fprintf(stderr,"Records up to %d already exist in log file\n",max_index);
    exit(1);
  }

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  /* Generate new states and write them to the log file. */

  n = dist_count_vars();
  q = chk_alloc (n, sizeof(double));

  for ( ; index<=max_index; index++)
  {
    dist_sample_prior (dst, q);

    logf.header.type = 'q';
    logf.header.index = index;
    logf.header.size = n * sizeof (double);
    log_file_append (&logf, q);

    logf.header.type = 'r';
    logf.header.index = index;
    logf.header.size = sizeof (rand_state);
    log_file_append (&logf, rand_get_state());
  }

  log_file_close(&logf);

  exit(0);
}
