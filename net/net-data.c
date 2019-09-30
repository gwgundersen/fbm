/* NET-DATA.C - Procedures for reading data for neural networks. */

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
#include "net.h"
#include "net-data.h"
#include "numin.h"


/* VARIABLES HOLDING DATA.  As declared in net-data.h. */

data_specifications *data_spec;	/* Specifications of data sets */

int N_train;			/* Number of training cases */

net_values *train_values;	/* Values associated with training cases */
double *train_targets;		/* True targets for training cases */

int N_test;			/* Number of test cases */

net_values *test_values;	/* Values associated with test cases */
double *test_targets;		/* True targets for test cases */


/* PROCEDURES. */

static double     *read_targets (numin_source *, int,   net_arch *);
static net_values *read_inputs  (numin_source *, int *, net_arch *, 
                                 model_specification *, model_survival *);


/* FREE SPACE OCCUPIED BY DATA.  Also useful as a way to reset when the
   architecture changes. */

void net_data_free (void)
{
  if (train_values!=0)
  { free(train_values);
    train_values = 0;
    N_train = 0;
  }

  if (train_targets!=0)
  { free(train_targets);
    train_targets = 0;
  }

  if (test_values!=0)
  { free(test_values);
    test_values = 0;
    N_test = 0;
  }

  if (test_targets!=0)
  { free(test_targets);
    test_targets = 0;
  }
}


/* READ TRAINING AND/OR TEST DATA FOR NETWORK.  Either or both of these
   data sets are read, depending on the options passed.  If a data set
   has already been read, it isn't read again.  This procedure also checks
   that the data specifications are consistent with the network architecture. 

   For survival models with non-constant hazard, the first input in a case, 
   representing time, is set to zero by this procedure. */

void net_data_read
( int want_train,	/* Do we want the training data? */
  int want_test,	/* Do we want the test data? */
  net_arch *arch,	/* Network architecture */
  model_specification *model, /* Data model being used */
  model_survival *surv	/* Survival model, or zero if irrelevant */
)
{
  numin_source ns;

  if (train_values!=0) want_train = 0;
  if (test_values!=0)  want_test = 0;

  if (model_targets(model,arch->N_outputs) != data_spec->N_targets
   || arch->N_inputs != data_spec->N_inputs 
                    + (model!=0 && model->type=='V' && surv->hazard_type!='C'))
  { fprintf(stderr,
     "Number of inputs/targets in data specification doesn't match network\n");
    exit(1);
  }

  if (model!=0 && model->type=='C' && arch->N_outputs!=data_spec->int_target)
  { fprintf(stderr,
"Integer range for targets does not match number of outputs for class model\n");
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
    train_values = read_inputs (&ns, &N_train, arch, model, surv);

    numin_spec (&ns, data_spec->train_targets, data_spec->N_targets);
    train_targets = read_targets (&ns, N_train, arch);
  }

  if (want_test && data_spec->test_inputs[0]!=0)
  {
    numin_spec (&ns, "data@1,0",1);
    numin_spec (&ns, data_spec->test_inputs, data_spec->N_inputs);
    test_values = read_inputs (&ns, &N_test, arch, model, surv);

    if (data_spec->test_targets[0]!=0)
    { numin_spec (&ns, data_spec->test_targets, data_spec->N_targets);
      test_targets = read_targets (&ns, N_test, arch);
    }
  }
}


/* READ INPUTS VALUES FOR A SET OF CASES. */

static net_values *read_inputs
( numin_source *ns,
  int *N_cases_ptr,
  net_arch *arch,
  model_specification *model, 
  model_survival *surv
)
{
  net_value *value_block;
  net_values *values;
  int value_count;
  int N_cases;
  int i, j, j0;

  N_cases = numin_start(ns);

  value_count = net_setup_value_count(arch);

  value_block = chk_alloc (value_count*N_cases, sizeof *value_block);
  values      = chk_alloc (N_cases, sizeof *values);

  for (i = 0; i<N_cases; i++) 
  { net_setup_value_pointers (&values[i], value_block+value_count*i, arch);
  }

  for (i = 0; i<N_cases; i++) 
  { if (model!=0 && model->type=='V' && surv->hazard_type!='C')
    { values[i].i[0] = 0;
      j0 = 1;
    }
    else
    { j0 = 0;
    }
    numin_read(ns,values[i].i+j0);
    for (j = j0; j<arch->N_inputs; j++)
    { values[i].i[j] 
        = data_trans (values[i].i[j], data_spec->trans[j-j0]);
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
  net_arch *arch
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
