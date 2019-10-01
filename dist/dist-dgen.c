/* DIST-DGEN.C - Program to sample targets given parameter values. */

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
#include "dist.h"
#include "formula.h"
#include "rand.h"
#include "data.h"
#include "numin.h"
#include "dist-data.h"

static void usage (void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  char *targets_file;
  char *inputs;
  int index;

  char **ip;
  int c, i, j;

  log_file logf;
  log_gobbled logg;

  dist_spec *dst;

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

  if (0)
  { printf("%d %s %s\n",index,inputs,targets_file);
  }

  /* Open log file and read specification. */

  log_file_open (&logf, index<0);

  log_gobble_init(&logg,0);
  logg.req_size['d'] = sizeof *dst;
  logg.req_size['r'] = sizeof (rand_state);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  /* Look at specification, to find out what state variables exist. */

  if ((dst = logg.data['d'])==0)
  { fprintf(stderr,"No distribution specification in log file\n");
    exit(1);
  }

  if (!dst->Bayesian)
  { fprintf(stderr,
      "Generation of targets can be done only with Bayesian models\n");
    exit(1);
  }

  dist_const_eval(dst);  /* Also parses prior and likelihood formulas */

  if (formula_var_exists['t'-'a'][0] && formula_var_exists['t'-'a'][10])
  { fprintf(stderr,"Only one of the synonyms 't' and 't0' may be used\n");
    exit(1);
  }

  /* Read records in log file up to the desired iteration, and set random
     number state left after that iteration was generated. */

  while (!logf.at_end && (index<0 || logf.header.index<=index))
  { log_gobble(&logf,&logg);
  }

  if (logf.header.index<0 || index>=0 && logg.last_index!=index)
  { fprintf(stderr,"Requested index does not exist in log file\n");
    exit(1);
  }

  if (logg.data['q']==0 || logg.index['q'] != logg.last_index)
  { fprintf(stderr,"Missing parameter values for specified iteration\n");
    exit(1);
  }

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  /* Read inputs. */

  data_spec = logg.data['D'];

  if (data_spec==0) 
  { fprintf(stderr,"No data specification record\n");
    exit(1);
  }

  if (logg.actual_size['D'] != 
      data_spec_size(data_spec->N_inputs,data_spec->N_targets))
  { fprintf(stderr,"Data specification record is the wrong size!\n");
    exit(1);
  }

  if (data_spec->N_targets>10)
  { fprintf(stderr,"Can't have more then 10 target values\n");
    exit(1);
  }

  if (data_spec->N_inputs>10)
  { fprintf(stderr,"Can't have more then 10 input values\n");
    exit(1);
  }

  if (inputs!=0)
  { strcpy(data_spec->train_inputs,inputs);
  }

  numin_spec (&ns, "data@1,0",1);
  numin_spec (&ns, data_spec->train_inputs, data_spec->N_inputs);
  train_inputs = dist_data_read_inputs (&ns, &N_train);

  /* Set variables from state at specified iteration. */

  dist_unpack_vars (logg.data['q']);

  /* Open file to write targets to. */

  tf = fopen(targets_file,"w");
  if (tf==NULL)
  { fprintf(stderr,"Can't create targets file: %s\n",targets_file);
    exit(1);
  }

  /* Generate target values for each "training" case. */

  for (c = 0; c<N_train; c++)
  { 
    /* Set inputs to values for this case. */

    for (i = 0; i<data_spec->N_inputs; i++)
    { formula_var['i'-'a'][i] = train_inputs [data_spec->N_inputs*c + i];
    }
    formula_var['i'-'a'][10] = formula_var['i'-'a'][0];

    /* Sample values for the target variables. */

    formula_sample (dst->energy, "t");

    /* Consider 't' to be a synonym for 't0'. */

    if (formula_var_exists['t'-'a'][10])
    { formula_var['t'-'a'][0] = formula_var['t'-'a'][10];
    }

    /* Write out generated target values */

    for (j = 0; j<data_spec->N_targets; j++)
    { fprintf(tf," %16.10e",formula_var['t'-'a'][j]);
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
{ fprintf (stderr,
    "Usage: dist-dgen log-file [ index ] / [ inputs ] targets-file\n");
  exit(1);
}
