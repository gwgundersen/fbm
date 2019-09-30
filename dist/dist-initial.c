/* DIST-INITIAL.C - Set initial state. */

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
  dist_spec dspec, *dst = &dspec;
  mc_value *q;

  log_file logf;
  log_gobbled logg;

  int n_vars;
  int a, c, i;
  char *e, *p;

  /* Open log file and gobble up initial records. */

  if (argc<2) usage();
  logf.file_name = argv[1];
  log_file_open(&logf,1);

  log_gobble_init(&logg,0);
  logg.req_size['d'] = sizeof *dst;

  if (!logf.at_end && logf.header.index==-1)
  { log_gobble(&logf,&logg);
  }

  /* Look at specification, to find out what state variables exist. */

  if ((dst = logg.data['d'])==0)
  { fprintf(stderr,"No distribution specification in log file\n");
    exit(1);
  }

  (void) formula (dst->energy, 1, 0, 0);
  if (dst->Bayesian)
  { (void) formula (dst->energy + strlen(dst->energy) + 1, 1, 0, 0);
  }

  /* Make assignments to state variables as specified in the arguments.  The
     values are not allowed to be general formulas, since it's hard to 
     arrange for variable references to be prohibited. */

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
  }

  /* Put variables into MC state vector. */

  n_vars = dist_count_vars();
  q = chk_alloc (n_vars, sizeof(mc_value));
  dist_pack_vars(q);

  /* Append new state to log file. */

  logf.header.type = 'q';
  logf.header.index = 0;
  logf.header.size = n_vars * sizeof(mc_value);
  log_file_append(&logf,q);

  log_file_close(&logf);

  exit(0);

}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr, "Usage: dist-initial log-file { var=value }\n");

  exit(1);
}
