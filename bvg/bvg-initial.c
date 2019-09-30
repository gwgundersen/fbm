/* BVG-INITIAL.C - Set initial state for bivariate Gaussian. */

/* Copyright (c) 1996 by Radford M. Neal 
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
#include "log.h"
#include "mc.h"
#include "bvg.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
) 
{
  bvg_spec bspec, *bs = &bspec;
  mc_value *q;

  log_file logf;
  log_gobbled logg;

  int i;

  /* Open log file and gobble up initial records. */

  if (argc<4) usage();
  logf.file_name = argv[1];
  log_file_open(&logf,1);

  log_gobble_init(&logg,0);
  logg.req_size['B'] = sizeof *bs;

  if (!logf.at_end && logf.header.index==-1)
  { log_gobble(&logf,&logg);
  }

  if ((bs = logg.data['B'])==0)
  { fprintf(stderr,"No bivariate Gaussian specification in log file\n");
    exit(1);
  }

  /* Allocate space for state. */

  q = chk_alloc (2*bs->rep, sizeof(mc_value));

  /* Set state from arguments. */

  if (argc!=2+2*bs->rep) usage();

  for (i = 0; i<2*bs->rep; i++)
  { q[i] = atof(argv[2+i]);
  }

  /* Append new state to log file. */

  logf.header.type = 'X';
  logf.header.index = 0;
  logf.header.size = 2 * bs->rep * sizeof(mc_value);
  log_file_append(&logf,q);

  log_file_close(&logf);

  exit(0);

}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
    "Usage: bvg-initial log-file value-1 value-2 { value-1 value-2 }\n");

  exit(1);
}
