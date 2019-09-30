/* NET-GEN.C - Program to generate networks (eg, from prior distribution). */

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
#include "net.h"
#include "rand.h"


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  net_arch *a;
  net_flags *flgs;
  net_priors *p;

  model_specification *m;
  model_survival *v;

  net_sigmas  sigmas, *s = &sigmas;
  net_params  params, *w = &params;

  int max_index, index, fix;
  double value, out_value;

  log_file logf;
  log_gobbled logg;
 
  /* Look at program arguments. */

  max_index = 0;
  value = out_value = 0;

  fix = argc>2 && strcmp(argv[2],"fix")==0 ? 2
      : argc>3 && strcmp(argv[3],"fix")==0 ? 3
      : argc;

  if (argc<2 || argc>fix+3
   || fix>2 && (max_index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0
   || fix+1<argc && (value = out_value = atof(argv[fix+1]))<=0
   || fix+2<argc && (out_value = atof(argv[fix+2]))<=0)
  { fprintf(stderr,
"Usage: net-gen log-file [ max-index ] [ \"fix\" [ value [ out-value ] ] ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read network architecture and priors. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  net_record_sizes(&logg);
  logg.req_size['r'] = sizeof (rand_state);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  a = logg.data['A'];
  m = logg.data['M'];
  p = logg.data['P'];
  v = logg.data['V'];

  flgs = logg.data['F'];

  net_check_specs_present(a,p,0,m,v);

  /* Allocate space for parameters and hyperparameters. */

  s->total_sigmas = net_setup_sigma_count(a,flgs,m);
  w->total_params = net_setup_param_count(a,flgs);

  s->sigma_block = chk_alloc (s->total_sigmas, sizeof (net_sigma));
  w->param_block = chk_alloc (w->total_params, sizeof (net_param));

  net_setup_sigma_pointers (s, a, flgs, m);
  net_setup_param_pointers (w, a, flgs);

  /* Read last records in log file to see where to start, and to get random
     number state left after last network was generated. */

  index = log_gobble_last(&logf,&logg);

  if (logg.last_index<0) 
  { index = 0;
  }

  if (index>max_index)
  { fprintf(stderr,"Networks up to %d already exist in log file\n",max_index);
    exit(1);
  }

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  /* Generate new networks and write them to the log file. */

  for ( ; index<=max_index; index++)
  {
    net_prior_generate (w, s, a, flgs, m, p, argv[fix]!=0, value, out_value);

    logf.header.type = 'S';
    logf.header.index = index;
    logf.header.size = s->total_sigmas * sizeof (net_sigma);
    log_file_append (&logf, s->sigma_block);

    logf.header.type = 'W';
    logf.header.index = index;
    logf.header.size = w->total_params * sizeof (net_param);
    log_file_append (&logf, w->param_block);

    if (argv[fix]==0)
    { logf.header.type = 'r';
      logf.header.index = index;
      logf.header.size = sizeof (rand_state);
      log_file_append (&logf, rand_get_state());
    }
  }

  log_file_close(&logf);

  exit(0);
}
