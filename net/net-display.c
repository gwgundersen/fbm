/* NET-DISPLAY.C - Program to print network parameters and other such info. */

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


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  net_arch *a;
  net_flags *flgs;
  model_specification *m;

  net_params  params, *w = &params;
  net_sigmas  sigmas, *s = &sigmas;

  log_file logf;
  log_gobbled logg;

  int sigmas_only;
  int params_only;
  int index;

  /* Look at arguments. */

  sigmas_only = 0;
  params_only = 0;

  if (argc>1 && (strcmp(argv[1],"-p")==0 || strcmp(argv[1],"-w")==0))
  { params_only = 1;
    argv += 1;
    argc -= 1;
  }
  else if (argc>1 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"-s")==0))
  { sigmas_only = 1;
    argv += 1;
    argc -= 1;
  }

  index = -1;

  if (argc!=2 && argc!=3 
   || argc>2 && (index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0) 
  { fprintf(stderr,"Usage: net-display [ -p | -h ] log-file [ index ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read network architecture. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  net_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  a = logg.data['A'];
  m = logg.data['M'];
  flgs = logg.data['F'];
  
  if (a==0)
  { fprintf(stderr,"No architecture specification in log file\n");
    exit(1);
  }

  s->total_sigmas = net_setup_sigma_count(a,flgs,m);
  w->total_params = net_setup_param_count(a,flgs);

  logg.req_size['S'] = s->total_sigmas * sizeof(net_sigma);
  logg.req_size['W'] = w->total_params * sizeof(net_param);

  /* Read the desired network from the log file. */

  if (index<0)
  { 
    log_gobble_last(&logf,&logg);

    if (logg.last_index<0)
    { fprintf(stderr,"No network in log file\n");
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
    { fprintf(stderr,"No network with that index is in the log file\n");
      exit(1);
    }

    log_gobble(&logf,&logg);
  }

  if (logg.index['W']!=index)
  { fprintf(stderr,"No weights stored for the network with that index\n");
    exit(1);
  }

  if (logg.index['S']!=index)
  { fprintf(stderr,"No sigmas stored for the network with that index\n");
    exit(1);
  }

  s->sigma_block = logg.data['S'];
  w->param_block = logg.data['W'];

  net_setup_sigma_pointers (s, a, flgs, m);
  net_setup_param_pointers (w, a, flgs);

  /* Print values of the parameters and hyperparameters, or whatever. */

  printf("\nNetwork in file \"%s\" with index %d%s\n", logf.file_name, index,
    sigmas_only ? " (sigmas only)" : params_only ? " (parameters only)" : "");

  if (params_only)
  { net_print_params(w,0,a,flgs,m);
  }
  else if (sigmas_only)
  { net_print_sigmas(s,a,flgs,m);
  }
  else
  { net_print_params(w,s,a,flgs,m);
  }

  printf("\n");

  exit(0);
}
