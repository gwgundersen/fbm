/* MIX-DATA.C - Procedures for reading data for mixture model. */

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
#include "data.h"
#include "prior.h"
#include "model.h"
#include "numin.h"
#include "mix.h"
#include "mix-data.h"


/* VARIABLES HOLDING DATA.  As declared in mix-data.h. */

data_specifications *data_spec;	/* Specifications of data sets */

int N_train;			/* Number of training cases */

double *train_inputs;		/* Inputs for training cases */
double *train_targets;		/* True targets for training cases */

int N_test;			/* Number of test cases */

double *test_inputs;		/* Inputs for test cases */
double *test_targets;		/* True targets for test cases */


/* PROCEDURES. */

static double *read_targets (numin_source *, int, mix_spec *);
static double *read_inputs  (numin_source *, int *, mix_spec *);


/* FREE SPACE OCCUPIED BY DATA.  Also useful as a way to reset when the
   mixture model specification changes. */

void mix_data_free (void)
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
   specifications are consistent with the mixture model specifications. */

void mix_data_read
( int want_train,	/* Do we want the training data? */
  int want_test,	/* Do we want the test data? */
  mix_spec *mx,		/* Mixture model specification */
  model_specification *model /* Data model being used */
)
{
  numin_source ns;

  if (train_inputs!=0) want_train = 0;
  if (test_inputs!=0)  want_test = 0;

  if (mx->N_targets!=data_spec->N_targets 
   || mx->N_inputs != data_spec->N_inputs)
  { fprintf(stderr,
     "Number of inputs/targets in data specification doesn't match model\n");
    exit(1);
  }

  if (model!=0 && model->type=='B' && data_spec->int_target!=2)
  { fprintf(stderr,"Data for binary targets must be specified to be binary\n");
    exit(1);
  }

  if (want_train)
  { 
    numin_spec (&ns, "data@1,0",1);
    numin_spec (&ns, data_spec->train_inputs, data_spec->N_inputs);
    train_inputs = read_inputs (&ns, &N_train, mx);

    numin_spec (&ns, data_spec->train_targets, data_spec->N_targets);
    train_targets = read_targets (&ns, N_train, mx);
  }

  if (want_test && data_spec->test_inputs[0]!=0)
  {
    numin_spec (&ns, "data@1,0",1);
    numin_spec (&ns, data_spec->test_inputs, data_spec->N_inputs);
    test_inputs = read_inputs (&ns, &N_test, mx);

    if (data_spec->test_targets[0]!=0)
    { numin_spec (&ns, data_spec->test_targets, data_spec->N_targets);
      test_targets = read_targets (&ns, N_test, mx);
    }
  }
}


/* READ INPUTS VALUES FOR A SET OF CASES. */

static double *read_inputs
( numin_source *ns,
  int *N_cases_ptr,
  mix_spec *mx
)
{
  double *values;
  int N_cases;
  int i, j;

  N_cases = numin_start(ns);

  values = chk_alloc (data_spec->N_inputs*N_cases, sizeof (double));

  for (i = 0; i<N_cases; i++) 
  { numin_read(ns,values+data_spec->N_inputs*i);
    for (j = 0; j<mx->N_inputs; j++)
    { values[data_spec->N_inputs*i+j] =
       data_trans (values[data_spec->N_inputs*i+j] , data_spec->trans[j]);
    }
  }

  numin_close(ns);

  *N_cases_ptr = N_cases;

  return values;
}


/* READ TARGET VALUES FOR A SET OF CASES. */

static double *read_targets
( numin_source *ns,
  int N_cases,
  mix_spec *mx
)
{
  double *tg;
  int i, j;

  if (numin_start(ns)!=N_cases)
  { fprintf(stderr,
      "Number of input cases doesn't match number of target cases\n");
    exit(1);
  }

  tg = chk_alloc (data_spec->N_targets*N_cases, sizeof (double));

  for (i = 0; i<N_cases; i++)
  { 
    numin_read(ns,tg+data_spec->N_targets*i);

    for (j = 0; j<data_spec->N_targets; j++)
    { tg[data_spec->N_targets*i+j] =
         data_trans (tg[data_spec->N_targets*i+j], 
                     data_spec->trans[data_spec->N_inputs+j]);
    }
  }

  numin_close(ns);

  return tg;
}
