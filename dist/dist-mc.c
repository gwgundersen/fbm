/* DIST-MC.C - Markov chain Monte Carlo for a specified distribution. */

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
#include "formula.h"
#include "data.h"
#include "dist.h"
#include "dist-data.h"


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* LOCAL VARIABLES. */

static dist_spec *dst;	/* Specification of distribution */
static double *ss;	/* Stepsizes selected, or zero. */
static int uses_data;	/* Does likelihood refer to data vars ("i" or "t")? */


/* SET UP REQUIRED RECORD SIZES PRIOR TO GOBBLING RECORDS. */

void mc_app_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{ logg->req_size['d'] = sizeof (dist_spec);
}


/* INITIALIZE AND SET UP DYNAMIC STATE STRUCTURE.  Also evaluates any
   constants used in the distribution specification. */

void mc_app_initialize
( log_gobbled *logg,	/* Records gobbled up from head and tail of log file */
  mc_dynamic_state *ds	/* Structure holding pointers to dynamical state */
)
{ 
  char *a, *f, *p;
  int c, i;

  if ((dst = logg->data['d'])==0)
  { fprintf(stderr,"No distribution specification in log file\n");
    exit(1);
  }

  a = dst->energy + strlen(dst->energy) + 1;
  if (dst->Bayesian) a += strlen(a) + 1;

  while (*a)
  { f = formula_def(a,&c,&i);
    formula_var[c][i] = formula(f,0,1,0);
    formula_var_exists[c][i] = 1;
    a += strlen(a) + 1;
  }

  (void) formula (dst->energy, 1, 0, 0);

  if (dst->Bayesian)
  { 
    (void) formula (dst->energy + strlen(dst->energy) + 1, 1, 0, 0);

    uses_data = 0;

    for (p = "it"; *p; p++)
    { for (i = 0; i<=10; i++)
      { if (formula_var_exists[*p-'a'][i])
        { uses_data = 1;
          break;
        }
      }
    }
  }

  ds->dim = dist_count_vars();

  ds->temp_state = 0;

  logg->req_size['q'] = ds->dim * sizeof(mc_value);

  if ((ds->q = logg->data['q'])==0)
  { 
    ds->q = chk_alloc (ds->dim, sizeof(mc_value));
    for (i = 0; i<ds->dim; i++) ds->q[i] = 0;
  }
  else
  {
    if (logg->index['q'] != logg->last_index)
    { fprintf(stderr,"Missing current point in iteration stored in log file\n");
      exit(1);
    }
  }

  ds->stepsize = chk_alloc (ds->dim, sizeof(mc_value));

  ss = logg->data['s'];

  if (ss!=0)
  { if (logg->actual_size['s']!=ds->dim*sizeof(double))
    { fprintf(stderr,"Stepsize record is the wrong size!\n");
      exit(1);
    }
    for (i = 0; i<ds->dim; i++) 
    { if (ss[i]<=0)
      { fprintf(stderr,"Stepsize record is garbled!\n");
        exit(1);
      }
    }
  }

  mc_app_stepsizes(ds);

  data_spec = logg->data['D'];

  if (data_spec && logg->actual_size['D'] !=
                     data_spec_size(data_spec->N_inputs,data_spec->N_targets))
  { fprintf(stderr,"Data specification record is the wrong size!\n");
    exit(1);
  }

  if (data_spec!=0)
  { 
    if (!dst->Bayesian || !uses_data)
    { 
      fprintf (stderr,
       "Warning: A data specification is inappropriate, and will be ignored\n");
    }
    else
    {
      dist_data_free();   
      dist_data_read();
    }
  }
}


/* SAVE POSITION AND AUXILIARY PART OF STATE. */

void mc_app_save
( mc_dynamic_state *ds,	/* Current dyanamical state */
  log_file *logf,	/* Log file state structure */
  int index		/* Index of iteration being saved */
)
{
  logf->header.type = 'q';
  logf->header.index = index;
  logf->header.size = ds->dim * sizeof(mc_value);
  log_file_append(logf,ds->q);
}


/* APPLICATION-SPECIFIC SAMPLING PROCEDURE. */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double a,
  double a2,
  mc_iter *it,
  mc_temp_sched *sch
)
{
  return 0;
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
  double sumsq, inv_temp, e;
  double *g;
  char *p;
  int i;

  inv_temp = !ds->temp_state ? 1 : ds->temp_state->inv_temp;

  if (dst->Bayesian) /* Bayesian model */
  {
    dist_unpack_vars(ds->q);

    if (inv_temp>=0)
    { e = dist_prior (dst, grad);
    }
    else
    { e = 0;
      if (grad)
      { for (i = 0; i<ds->dim; i++) 
        { grad[i] = 0;
        }
      }
      inv_temp = -inv_temp;
    }

    if (!uses_data) /* Bayesian model with data in likelihood formula itself */
    { 
      if (inv_temp!=0) 
      {
        e += inv_temp * formula (dst->energy, 0, 1, grad ? STATE_VARS : 0);

        if (grad)
        {
          g = grad;

          for (p = STATE_VARS; *p; p++)
          { for (i = 0; i<=10; i++)
            { if (formula_var_exists[*p-'a'][i])
              { *g++ += inv_temp * formula_gradient[*p-'a'][i]; 
              }
            }
          }
        }
      }
    }

    else if (N_train>0) /* Bayesian model with data in file */
    { 
      if (inv_temp!=0) e += dist_total_likelihood (dst, inv_temp, grad);
    }

    if (e>1e100)
    { e = 1e100;
      if (grad)
      { for (i = 0; i<ds->dim; i++)
        { grad[i] = 0;
        }
      }
    }

    if (energy) *energy = e;
  }

  else /* Not a Bayesian model */
  {
    if (inv_temp<0)
    { fprintf(stderr,
      "Hamiltonian importance sampling is possible only for Bayesian models\n");
      exit(1);
    }

    if (inv_temp!=1)
    { sumsq = 0;
      for (i = 0; i<ds->dim; i++)
      { sumsq += ds->q[i]*ds->q[i];
      }
    }
  
    if (inv_temp==0 && !dst->zero_temper)
    { if (energy) 
      { *energy = sumsq/2 + ds->dim*log(2*M_PI)/2;
      }
      if (grad)
      { for (i = 0; i<ds->dim; i++)
        { grad[i] = ds->q[i];
        }
      }
      return;
    }
  
    dist_unpack_vars(ds->q);
    e = formula (dst->energy, 0, 1, grad ? STATE_VARS : 0);
  
    if (energy) 
    { *energy = e>1e100 ? 1e100 : e;
      if (inv_temp!=1)
      { if (dst->zero_temper)
        { *energy *= inv_temp; 
        }
        else
        { *energy = *energy * inv_temp
            + (sumsq/2 + ds->dim*log(2*M_PI)/2) * (1-inv_temp);
        }
      }
    }
    
    if (grad)
    { if (e>1e100)
      { for (i = 0; i<ds->dim; i++) 
        { grad[i] = 0;
        }
      }
      else
      { dist_pack_grad(grad);
        if (inv_temp!=1)
        { if (dst->zero_temper)
          { for (i = 0; i<ds->dim; i++) 
            { grad[i] *= inv_temp;
            }
          } 
          else
          { for (i = 0; i<ds->dim; i++) 
            { grad[i] = grad[i] * inv_temp + ds->q[i] * (1-inv_temp);
            }
          } 
        }
      }
    }
  }
}


/* SAMPLE FROM DISTRIBUTION AT INVERSE TEMPERATURE OF ZERO.  Returns zero
   if this is not possible. */

int mc_app_zero_gen
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
  int i;

  /* For a Bayesian model, use the dist_sample_prior procedure. */

  if (dst->Bayesian)
  { 
    dist_sample_prior(dst,ds->q);

    return 1;
  }

  /* When -zero-temp was specified, sampling isn't possible. */

  else if (dst->zero_temper)
  { 
    return 0;
  }

  /* Otherwise, we sample from a standard Gaussian. */

  else
  { 
    for (i = 0; i<ds->dim; i++)
    { ds->q[i] = rand_gaussian();
    }

    return 1;
  }

}


/* SET STEPSIZES FOR EACH COORDINATE. */

void mc_app_stepsizes
( mc_dynamic_state *ds	/* Current dyanamical state */
)
{ 
  int i;

  if (ss==0)
  { for (i = 0; i<ds->dim; i++) 
    { ds->stepsize[i] = 1.0;
    }
  }
  else
  { for (i = 0; i<ds->dim; i++) 
    { ds->stepsize[i] = ss[i];
    }
  }
}
