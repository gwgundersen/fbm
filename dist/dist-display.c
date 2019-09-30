/* DIST-DISPLAY.C - Display state variables. */

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


static void usage(void)
{ fprintf (stderr, "Usage: dist-display log-file [ index ]\n");
  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  dist_spec *dst;
  log_file logf;
  log_gobbled logg;
  mc_value *q;
  int index;

  /* Look at arguments. */

  index = -1;

  if (argc!=2 && argc!=3 
   || argc>2 && (index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0) 
  { usage();
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  logg.req_size['d'] = sizeof *dst;

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  dst = logg.data['d'];
  
  if (dst==0)
  { fprintf(stderr,"No distribution specification in log file\n");
    exit(1);
  }

  (void) formula (dst->energy, 1, 0, 0);
  if (dst->Bayesian)
  { (void) formula (dst->energy + strlen(dst->energy) + 1, 1, 0, 0);
  }

  logg.req_size['q'] = dist_count_vars() * sizeof(mc_value);

  /* Read the desired records from the log file. */

  if (index<0)
  { 
    log_gobble_last(&logf,&logg);

    if (logg.last_index<0)
    { fprintf(stderr,"No state variable record in log file\n");
      exit(1);
    }

    index = logg.last_index;
  }
  else
  {
    while (!logf.at_end && logf.header.index!=index)
    { log_file_forward(&logf);
    }

    if (logf.at_end)
    { fprintf(stderr,"That index does not appear in the log file\n");
      exit(1);
    }

    log_gobble(&logf,&logg);
  }

  if (logg.index['q']!=index)
  { fprintf(stderr,"No state variables stored with that index\n");
    exit(1);
  }
  
  /* Unpack the state variables. */

  q = logg.data['q'];
  dist_unpack_vars(q);

  /* Print values of all state variables. */

  printf("\nSTATE VARIABLES IN FILE \"%s\" WITH INDEX %d\n\n", 
         logf.file_name, index);

  dist_print_vars();

  printf("\n");

  exit(0);
}
