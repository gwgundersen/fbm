/* BVG-MC.C - Bivariate Gaussian module for Markov chain Monte Carlo. */

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
#include "rand.h"
#include "log.h"
#include "mc.h"
#include "bvg.h"


/* LOCAL VARIABLES. */

static bvg_spec *bs;	/* Specification of bivariate Gaussian distribution */


/* SET UP REQUIRED RECORD SIZES PRIOR TO GOBBLING RECORDS. */

void mc_app_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{ logg->req_size['B'] = sizeof *bs;
}


/* INITIALIZE AND SET UP DYNAMIC STATE STRUCTURE. */

void mc_app_initialize
( log_gobbled *logg,	/* Records gobbled up from head and tail of log file */
  mc_dynamic_state *ds	/* Structure holding pointers to dynamical state */
)
{ 
  int i;

  if ((bs = logg->data['B'])==0)
  { fprintf(stderr,"No specification for bivariate Gaussian in log file\n");
    exit(1);
  }

  ds->dim = 2*bs->rep;

  ds->temp_state = 0;

  logg->req_size['X'] = 2 * bs->rep * sizeof(mc_value);

  if ((ds->q = logg->data['X'])==0)
  { 
    ds->q = chk_alloc (2*bs->rep, sizeof(mc_value));
    for (i = 0; i<2*bs->rep; i++) ds->q[i] = 0;
  }
  else
  {
    if (logg->index['X'] != logg->last_index)
    { fprintf(stderr,"Missing current point in iteration stored in log file\n");
      exit(1);
    }
  }

  ds->stepsize = chk_alloc (2*bs->rep, sizeof(mc_value));
  for (i = 0; i<2*bs->rep; i++) ds->stepsize[i] = 1;
}


/* SAVE POSITION AND AUXILIARY PART OF STATE. */

void mc_app_save
( mc_dynamic_state *ds,	/* Current dyanamical state */
  log_file *logf,	/* Log file state structure */
  int index		/* Index of iteration being saved */
)
{
  logf->header.type = 'X';
  logf->header.index = index;
  logf->header.size = 2 * bs->rep * sizeof(mc_value);
  log_file_append(logf,ds->q);
}


/* APPLICATION-SPECIFIC SAMPLING PROCEDURE.  Implements the "gibbs" operation,
   which does a Gibbs sampling iteration, and the "gibbs0" and "gibbs1"
   operations, which do overrelaxed updates for the two coordinates.  
   Returns zero for unknown operations. */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double a,
  double a2,
  mc_iter *it,
  mc_temp_sched *sch
)
{
  double tf, t, m;
  double std1, std2;
  int i;

  tf = ds->temp_state ? ds->temp_state->inv_temp : 1.0;

  std1 = bs->std1 / sqrt(tf);
  std2 = bs->std2 / sqrt(tf);

  if (strcmp(op,"gibbs")==0)
  {
    t = sqrt (1 - bs->corr * bs->corr);

    for (i = 0; i<2*bs->rep; i += 2)
    { ds->q[i+0] = bs->corr*std1*ds->q[i+1] / std2 
                    + t*std1 * rand_gaussian();

      ds->q[i+1] = bs->corr*std2*ds->q[i+0] / std1 
                    + t*std2 * rand_gaussian();
    }
 
    ds->know_grad = 0;
    ds->know_pot  = 0;

    return 1;
  }

  else if (strcmp(op,"gibbs0")==0)
  {
    t = sqrt (1 - bs->corr * bs->corr);

    for (i = 0; i<2*bs->rep; i += 2)
    { m = bs->corr*std1*ds->q[i+1] / std2;
      ds->q[i+0] = m + a * (ds->q[i+0] - m)
                    + sqrt(1-a*a) * t*std1 * rand_gaussian();
    }

    ds->know_grad = 0;
    ds->know_pot  = 0;

    return 1;
  }

  else if (strcmp(op,"gibbs1")==0)
  {
    t = sqrt (1 - bs->corr * bs->corr);

    for (i = 0; i<2*bs->rep; i += 2)
    { m = bs->corr*std2*ds->q[i+0] / std1;
      ds->q[i+1] = m + a * (ds->q[i+1] - m) 
                     + sqrt(1-a*a) * t*std2 * rand_gaussian();
    }

    ds->know_grad = 0;
    ds->know_pot  = 0;

    return 1;
  }

  else
  { return 0;
  }
}


/* EVALUATE POTENTIAL ENERGY AND ITS GRADIENT. */

void mc_app_energy
( mc_dynamic_state *ds,	/* Current dyanamical state */
  int N_approx,		/* Number of gradient approximations in use */
  int w_approx,		/* Which approximation to use this time */
  double *energy,	/* Place to store energy, null if not required */
  mc_value *grad	/* Place to store gradient, null if not required */
)
{
  double tf, f, v1, v2, c, E;
  double std1, std2;
  int i;

  if (N_approx!=1) abort();

  tf = ds->temp_state ? ds->temp_state->inv_temp : 1.0;

  std1 = bs->std1 / sqrt(tf);
  std2 = bs->std2 / sqrt(tf);

  f = 1 / (1 - bs->corr*bs->corr);

  v1 = std1*std1;
  v2 = std2*std2;

  c = bs->corr / (std1*std2);

  E = 0;
  for (i = 0; i<2*bs->rep; i+=2)
  { E += (f/2) * ( (ds->q[i+0]*ds->q[i+0]) / v1 + (ds->q[i+1]*ds->q[i+1]) / v2 
                     - 2*c * (ds->q[i+0]*ds->q[i+1]) );
  }

  if (E>=1e30) 
  { E = 1e30;
    if (grad!=0) 
    { for (i = 0; i<2*bs->rep; i++) grad[i] = 0;
    }
  }
  else 
  { if (grad!=0) 
    { for (i = 0; i<2*bs->rep; i+=2)
      { grad[i+0] = f * (ds->q[i+0] / v1 - c * ds->q[i+1]);
        grad[i+1] = f * (ds->q[i+1] / v2 - c * ds->q[i+0]);
      }
    }
  }
 
  if (energy!=0) *energy = E;
}


/* SAMPLE FROM DISTRIBUTION AT INVERSE TEMPERATURE OF ZERO.  Returns zero
   if this is not possible. */

int mc_app_zero_gen
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
  return 0;
}


/* SET STEPSIZES FOR EACH COORDINATE. */

void mc_app_stepsizes
( mc_dynamic_state *ds	/* Current dyanamical state */
)
{ 
  int i;

  for (i = 0; i<2*bs->rep; i++)
  { ds->stepsize[i] = ds->temp_state ? 1 / sqrt(ds->temp_state->inv_temp) : 1;
  }
}
