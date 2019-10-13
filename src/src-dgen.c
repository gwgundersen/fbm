/* SRC-DGEN.C - Randomly generate detector measurements given src parameters. */

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
#include "rand.h"
#include "data.h"
#include "numin.h"
#include "src.h"
#include "src-data.h"

static void usage (void);


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  char *targets_file;
  char *inputs;
  int index;

  char **ip;
  int c, i, j;
  double v;

  log_file logf;
  log_gobbled logg;

  src_spec *src;
  det_spec *det;
  flow_spec *flow;

  src_params *params;

  numin_source ns;
  FILE *tf;

  /* Look at program arguments. */

  if (argc<4 || argc>6) usage();

  logf.file_name = argv[1];

  if (strcmp(argv[3],"/")==0)
  { if ((index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0)
    { usage(); 
    }
    if (argc<5) usage();
    ip = argv+4;
  }
  else if (strcmp(argv[2],"/")==0)
  { index = -1;
    if (argc>5) usage();
    ip = argv+3;
  }
  else
  { usage();
  }

  if (ip[1]==0)
  { inputs = 0;
    targets_file = ip[0];
  }
  else
  { inputs = ip[0];
    targets_file = ip[1];
  }

  if (0) /* for debugging */
  { printf("%d %s %s\n",index,inputs,targets_file);
  }

  /* Open log file and read specifications. */

  log_file_open (&logf, index<0);

  log_gobble_init(&logg,0);
  src_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  src = logg.data['S'];
  det = logg.data['T'];
  flow = logg.data['F'];

  src_check_specs_present(src,det,flow);

  data_spec = logg.data['D'];

  if (data_spec==0)
  { fprintf(stderr,"No data specifications in log file\n");
    exit(1);
  }

  if (logg.actual_size['D'] !=
       data_spec_size(data_spec->N_inputs,data_spec->N_targets))
  { fprintf(stderr,"Data specification record is the wrong size!\n");
    exit(1);
  }

  /* Read records in log file up to the desired iteration, and check that
     things are as they should be. */

  logg.req_size['q'] = src_params_length(src->highN);

  while (!logf.at_end && (index<0 || logf.header.index<=index))
  { log_gobble(&logf,&logg);
  }

  if (logf.header.index<0 || index>=0 && logg.last_index!=index)
  { fprintf(stderr,"Requested index does not exist in log file\n");
    exit(1);
  }

  if (logg.data['q']==0 || logg.index['q'] != logg.last_index)
  { fprintf(stderr,"Missing parameters for specified iteration\n");
    exit(1);
  }

  params = (src_params *) logg.data['q'];

  /* Set random number state left after the specified iteration was generated.*/

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  /* Read inputs. */

  numin_spec (&ns, "data@1,0",1);
  numin_spec (&ns, inputs ? inputs : data_spec->train_inputs, 
                   data_spec->N_inputs);
  train_inputs = src_data_read_inputs (&ns, &N_train);

  /* Open file to write targets to. */

  tf = fopen(targets_file,"w");
  if (tf==NULL)
  { fprintf(stderr,"Can't create targets file: %s\n",targets_file);
    exit(1);
  }

  /* Generate target values for each "training" case. */

  for (c = 0; c<N_train; c++)
  { 
    v = src_total (src, flow, params, train_inputs+4*c);

    if (params->inv_df==0)
    { v += exp(params->log_width) * rand_gaussian();
    }
    else
    { v += exp(params->log_width) * rand_gaussian() 
             * sqrt ((0.5/params->inv_df) / rand_gamma(0.5/params->inv_df));
    }

    fprintf(tf,"%15.5e\n",v);
  }

  fclose(tf);

  /* Append random state to file if using last index by default. */

  if (index<0)
  {
    logf.header.type = 'r';
    logf.header.index = logg.last_index;
    logf.header.size = sizeof (rand_state);
    log_file_append (&logf, rand_get_state());
  }

  log_file_close(&logf);

  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

static void usage (void)
{ fprintf (stderr, 
    "Usage: src-dgen log-file [ index ] / [ inputs ] targets-file\n");
  exit(1);
}
