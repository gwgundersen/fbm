/* SRC-DATA.C - Procedures for reading data for source model. */

/* Copyright (c) 1995-2007 by Radford M. Neal 
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
#include "numin.h"
#include "src.h"
#include "data.h"
#include "src-data.h"


/* VARIABLES HOLDING DATA.  As declared in src-data.h. */

data_specifications *data_spec;	/* Specifications of data sets */

int N_train;			/* Number of training cases */

double *train_inputs;		/* Inputs for training cases */
double *train_targets;		/* True targets for training cases */

int N_test;			/* Number of test cases */

double *test_inputs;		/* Inputs for test cases */
double *test_targets;		/* True targets for test cases */


/* FREE SPACE OCCUPIED BY DATA.  */

void src_data_free (void)
{
  if (train_inputs!=0)
  { free(train_inputs);
    train_inputs = 0;
    N_train = 0;
  }

  if (train_targets!=0)
  { free(train_targets);
    train_targets = 0;
  }

  if (test_inputs!=0)
  { free(test_inputs);
    test_inputs = 0;
    N_test = 0;
  }

  if (test_targets!=0)
  { free(test_targets);
    test_targets = 0;
  }
}


/* READ TRAINING AND/OR TEST DATA.  Either or both of these data sets are 
   read, depending on the options passed.  If a data set has already been 
   read, it isn't read again.  This procedure also checks that the data 
   specifications are consistent with the source model specifications. */

void src_data_read
( int want_train,	/* Do we want the training data? */
  int want_test		/* Do we want the test data? */
)
{
  numin_source ns;

  if (train_inputs!=0) want_train = 0;
  if (test_inputs!=0)  want_test = 0;

  if (data_spec->N_inputs>4)
  { fprintf(stderr,
"Number of inputs in data specification must be 4 or less for source models\n");
    exit(1);
  }

  if (data_spec->N_targets!=1)
  { fprintf(stderr,
     "Number of targets in data specification must be 1 for source models\n");
    exit(1);
  }

  if (want_train)
  { 
    numin_spec (&ns, "data@1,0",1);
    numin_spec (&ns, data_spec->train_inputs, data_spec->N_inputs);
    train_inputs = src_data_read_inputs (&ns, &N_train);

    if (data_spec->train_targets[0]!=0)
    { numin_spec (&ns, data_spec->train_targets, data_spec->N_targets);
      train_targets = src_data_read_targets (&ns, N_train);
    }
  }

  if (want_test && data_spec->test_inputs[0]!=0)
  {
    numin_spec (&ns, "data@1,0",1);
    numin_spec (&ns, data_spec->test_inputs, data_spec->N_inputs);
    test_inputs = src_data_read_inputs (&ns, &N_test);

    if (data_spec->test_targets[0]!=0)
    { numin_spec (&ns, data_spec->test_targets, data_spec->N_targets);
      test_targets = src_data_read_targets (&ns, N_test);
    }
  }
}


/* READ INPUT VALUES FOR A SET OF CASES. */

double *src_data_read_inputs
( numin_source *ns,
  int *N_cases_ptr
)
{
  double *values;
  int N_cases;
  int i, j;

  N_cases = numin_start(ns);

  values = chk_alloc (4*N_cases, sizeof (double));

  for (i = 0; i<N_cases; i++) 
  { numin_read(ns,values+4*i);
    for (j = 0; j<data_spec->N_inputs; j++)
    { values[4*i+j] = data_trans (values[4*i+j], data_spec->trans[j]);
    }
    for (; j<4; j++)
    { values[4*i+j] = 0;
    }
  }

  numin_close(ns);

  *N_cases_ptr = N_cases;

  return values;
}


/* READ TARGET VALUES FOR A SET OF CASES. */

double *src_data_read_targets
( numin_source *ns,
  int N_cases
)
{
  double *tg;
  int i, j;

  if (numin_start(ns)!=N_cases)
  { fprintf(stderr,
      "Number of input cases doesn't match number of target cases\n");
    exit(1);
  }

  tg = chk_alloc (N_cases, sizeof (double));

  if (data_spec->trans[data_spec->N_inputs].take_log)
  { fprintf (stderr,
             "Log transformation is not allowed for target measurements\n");
    exit(1);
  }

  for (i = 0; i<N_cases; i++)
  { numin_read(ns,tg+i);
    tg[i] = data_trans (tg[i], data_spec->trans[data_spec->N_inputs]);
  }

  numin_close(ns);

  return tg;
}
