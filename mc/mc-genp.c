/* MC-GENP.C - Skeleton program for generating random values for the momentum.*/

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
#include "rand.h"
#include "log.h"
#include "mc.h"


/* THE FOLLOWING IS NEEED BY THE PLOTTING ROUTINES. */

enum { PLT, TBL, HIST } program_type;


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  mc_dynamic_state ds;
  mc_ops *ops;

  log_file logf;
  log_gobbled logg;
  int index;
 
  int i;

  /* Look at program arguments. */

  if (argc!=2)
  { fprintf(stderr,"Usage: xxx-genp log-file\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and gobble all records.  Also sets index to the index 
     of the last record, or zero. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  mc_record_sizes(&logg);

  while (!logf.at_end)
  { log_gobble(&logf,&logg);
  }

  index = logg.last_index;
  if (index<0)
  { index = 0;
  }

  /* Exit if we don't need momentum variables. */

  ops = logg.data['o'];

  if (ops==0) 
  { exit(0);
  }

  for (i = 0; ; i++)
  { 
    if (i==Max_mc_ops || ops->op[i].type==0) 
    { exit(0);
    }

    if (strchr(MC_needs_p,ops->op[i].type)!=0) 
    { break;
    }
  }

  /* Initialize using records found. */

  mc_app_initialize(&logg,&ds);

  /* Generate random momentum variables. */

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  ds.p = chk_alloc (ds.dim, sizeof *ds.p);

  for (i = 0; i<ds.dim; i++)
  { ds.p[i] = rand_gaussian();
  }

  /* Write out momentum variables. */

  logf.header.type = 'p';
  logf.header.index = index;
  logf.header.size = ds.dim * sizeof *ds.p;
  log_file_append (&logf, ds.p);

  /* Write out new random number generator state. */

  logf.header.type = 'r';
  logf.header.index = index;
  logf.header.size = sizeof (rand_state);
  log_file_append (&logf, rand_get_state());

  exit(0);
}
