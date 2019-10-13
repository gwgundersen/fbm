/* GP-DGEN.C - Generate values for target variables given latent values. */

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
#include "rand.h"
#include "data.h"
#include "model.h"
#include "gp.h"
#include "gp-data.h"

static void usage (void);


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  char *targets_file;
  int index;

  int c, i, j;
  double v;

  log_file logf;
  log_gobbled logg;

  gp_spec *gp;
  model_specification *m;
  data_specifications *data_spec;
  gp_hypers hypers, *h = &hypers;

  double *values, *variances, *nvar, *targets;

  int N_train;
  FILE *tf;

  /* Look at program arguments. */

  if (argc<4 || argc>5) usage();

  logf.file_name = argv[1];

  if (strcmp(argv[3],"/")==0)
  { if ((index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0)
    { usage(); 
    }
    if (argc<5) usage();
    targets_file = argv[4];
  }
  else if (strcmp(argv[2],"/")==0)
  { index = -1;
    if (argc>4) usage();
    targets_file = argv[3];
  }
  else
  { usage();
  }

  /* Open log file and read specifications. */

  log_file_open (&logf, index<0);

  log_gobble_init(&logg,0);
  gp_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  gp = logg.data['P'];
  m = logg.data['M'];
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

  gp_check_specs_present (gp,0,m,0);

  if (m!=0 && m->type=='V')
  { fprintf(stderr,"Can't handle survival models\n");
    exit(1);
  }

  model_values_check(m,data_spec,gp->N_outputs,0);

  /* Read records in log file up to the desired iteration, and check that
     things are as they should be. */

  h->total_hypers = gp_hyper_count(gp,m);
  logg.req_size['S'] = h->total_hypers * sizeof(double);

  while (!logf.at_end && (index<0 || logf.header.index<=index))
  { log_gobble(&logf,&logg);
  }

  if (logf.header.index<0 || index>=0 && logg.last_index!=index)
  { fprintf(stderr,"Requested index does not exist in log file\n");
    exit(1);
  }

  if (logg.data['S']==0 || logg.index['S'] != logg.last_index)
  { fprintf(stderr,"Missing hyperparameters for specified iteration\n");
    exit(1);
  }

  h->hyper_block = logg.data['S'];
  gp_hyper_pointers (h, gp, m);

  if (logg.data['F']==0 || logg.index['F'] != logg.last_index)
  { fprintf(stderr,"Missing latent values for specified iteration\n");
    exit(1);
  }

  values = (double *) logg.data['F'];
  N_train = logg.actual_size['F'] / sizeof(double) / gp->N_outputs;

  variances = logg.index['N']==logg.last_index ? logg.data['N'] : 0;

  if (variances!=0)
  { if (logg.actual_size['N'] != N_train*sizeof(double))
    { fprintf(stderr,"Record with noise variances is the wrong size!\n");
      exit(1);
    }
  }

  /* Set up noise variances on a case-by-case or fixed basis. */

  if (m!=0 && m->type=='R')
  { if (m->noise.alpha[2]!=0)
    { if (variances==0)
      { fprintf(stderr,
          "Missing case-by-case noise variances at specified iteration\n");
        exit(1);
      }
    }
    else 
    { variances = 0;
      nvar = chk_alloc(data_spec->N_targets,sizeof(double));
      for (j = 0; j<data_spec->N_targets; j++)
      { nvar[j] = exp (2 * *h->noise[j]);
      }
    }
  }
  else
  { variances = 0;
    nvar = 0;
  }

  /* Set random number state left after the specified iteration was generated.*/

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  /* Open file to write targets to. */

  tf = fopen(targets_file,"w");
  if (tf==NULL)
  { fprintf(stderr,"Can't create targets file: %s\n",targets_file);
    exit(1);
  }

  /* Generate target values for each "training" case. */

  targets = chk_alloc(data_spec->N_targets,sizeof(double));

  for (c = 0; c<N_train; c++)
  { 
    model_gen (m, data_spec, 
               variances ? variances+c*data_spec->N_targets : nvar,
               values+c*gp->N_outputs, gp->N_outputs,
               targets);

    for (j = 0; j<data_spec->N_targets; j++)
    { fprintf(tf," %16.10e",targets[j]);
    }

    fprintf(tf,"\n");
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
{ fprintf (stderr, "Usage: gp-dgen log-file [ index ] / targets-file\n");
  exit(1);
}
