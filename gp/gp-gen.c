/* GP-GEN.C - Program to generate Gaussian process hyperparameters. */

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
#include "model.h"
#include "data.h"
#include "gp.h"
#include "rand.h"


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  gp_spec *gp;

  model_specification *m;
  model_survival *v;

  gp_hypers hypers, *h = &hypers;

  int max_index, index, fixarg, fix;
  double scale_value, relevance_value;

  log_file logf;
  log_gobbled logg;
 
  /* Look at program arguments. */

  max_index = 0;
  scale_value = 0;

  fixarg = argc>2 && strcmp(argv[2],"fix")==0 ? 2
         : argc>3 && strcmp(argv[3],"fix")==0 ? 3
         : argc;

  fix = fixarg!=argc;

  if (argc<2 || argc!=fixarg && argc!=fixarg+1 && argc!=fixarg+3
   || fixarg>2 && (max_index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0
   || fixarg+1<argc && (scale_value = atof(argv[fixarg+1]))<=0
   || fixarg+2<argc && (relevance_value = atof(argv[fixarg+2]))<=0)
  { fprintf(stderr,
"Usage: gp-gen log-file [ max-index ] [ \"fix\" [ scale-value relevance-value ] ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  gp_record_sizes(&logg);
  logg.req_size['r'] = sizeof (rand_state);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  gp = logg.data['P'];
  m = logg.data['M'];
  v = logg.data['V'];

  gp_check_specs_present(gp,0,m,v);

  /* Allocate space for hyperparameters. */

  h->total_hypers = gp_hyper_count(gp,m);
  h->hyper_block = chk_alloc (h->total_hypers, sizeof (double));

  gp_hyper_pointers (h, gp, m);

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

  /* Generate new records with each index and write them to the log file. */

  for ( ; index<=max_index; index++)
  {
    /* Generate hyperparameter values. */

    gp_prior_generate (h, gp, m, fix, scale_value, relevance_value);

    /* Write out records, including random number state record if we've
       been using random numbers. */

    logf.header.type = 'S';
    logf.header.index = index;
    logf.header.size = h->total_hypers * sizeof (double);
    log_file_append (&logf, h->hyper_block);

    if (!fix)
    { logf.header.type = 'r';
      logf.header.index = index;
      logf.header.size = sizeof (rand_state);
      log_file_append (&logf, rand_get_state());
    }
  }

  log_file_close(&logf);

  exit(0);
}
