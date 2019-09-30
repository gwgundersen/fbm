/* DIST-DATA.C - Procedures for reading data for Bayesian model. */

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
#include "numin.h"
#include "dist-data.h"


/* VARIABLES HOLDING DATA.  As declared in dist-data.h. */

data_specifications *data_spec;	/* Specifications of data sets */

int N_train;			/* Number of training cases */

double *train_inputs;		/* Inputs for training cases */
double *train_targets;		/* True targets for training cases */


/* PROCEDURES. */

static double *read_inputs  (numin_source *, int *);
static double *read_targets (numin_source *, int);


/* FREE SPACE OCCUPIED BY DATA. */

void dist_data_free (void)
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
}


/* READ TRAINING DATA. If the data has already been read, it isn't read again.*/

void dist_data_read (void)
{
  numin_source ns;

  if (train_inputs!=0) return;

  numin_spec (&ns, "data@1,0",1);
  numin_spec (&ns, data_spec->train_inputs, data_spec->N_inputs);
  train_inputs = read_inputs (&ns, &N_train);

  numin_spec (&ns, data_spec->train_targets, data_spec->N_targets);
  train_targets = read_targets (&ns, N_train);
}


/* READ INPUTS VALUES FOR A SET OF CASES. */

static double *read_inputs
( numin_source *ns,
  int *N_cases_ptr
)
{
  double *values;
  int N_cases;
  int i, j;

  N_cases = numin_start(ns);

  values = chk_alloc (data_spec->N_inputs*N_cases, sizeof (double));

  for (i = 0; i<N_cases; i++) 
  { numin_read(ns,values+data_spec->N_inputs*i);
    for (j = 0; j<data_spec->N_inputs; j++)
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
