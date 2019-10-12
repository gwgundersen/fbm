/* SRC-GEN.C - Program to generate from prior for source location model. */

/* Copyright (c) 2007 by Radford M. Neal 
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
#include "src.h"
#include "rand.h"


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  log_file logf;
  log_gobbled logg;

  src_spec *src;
  det_spec *det;
  flow_spec *flow;

  src_params *params; 

  int index, max_index;
  char junk;
  int i, j;
 
  /* Look at program arguments. */

  max_index = 0;

  if (argc<2 || argc>3 || argc==3 && sscanf(argv[2],"%d%c",&max_index,&junk)!=1)
  { fprintf(stderr, "Usage: mix-gen log-file [ max-index ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  src_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  src = logg.data['S'];
  det = logg.data['T'];
  flow = logg.data['F'];

  src_check_specs_present (src, det, flow);

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

  /* Allocate space for parameters.  Allow enough space for the specified 
     maximum number of sources. */

  params = 
    chk_alloc (src_params_length(src->highN)/sizeof(double), sizeof(double));

  /* Generate new records with each index and write them to the log file. */

  for ( ; index<=max_index; index++)
  {
    /* Generate parameter values from prior. */

    src_prior_generate (params, src, det, flow);

    /* Write out parameter record, and random number state record. */

    logf.header.type = 'q';
    logf.header.index = index;
    logf.header.size = src_params_length(src->highN);
    log_file_append (&logf, params);

    logf.header.type = 'r';
    logf.header.index = index;
    logf.header.size = sizeof (rand_state);
    log_file_append (&logf, rand_get_state());
  }

  log_file_close(&logf);

  exit(0);
}
