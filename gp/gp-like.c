/* GP-LIKE.C - Likelihood computations for Gaussian Process models. */  
 
/* Copyright (c) 1998 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
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
#include "gp.h"
#include "rand.h"


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* COMPUTE LOG LIKELIHOOD FOR TARGETS GIVEN LATENT VALUES FOR ONE CASE. */

double gp_likelihood
( gp_hypers *h,			/* Hyperparameter values */
  model_specification *m,	/* Specification for model used */
  data_specifications *d,	/* Specification for data */
  double *targets,		/* Values of targets for this case */
  double *latent_values,	/* Values of latent variables for this case */
  double *noise_variances	/* Noise variances (if they vary by case) */
)
{ 
  double l, n, a;
  int i, j;
   
  if (m->type=='C')
  { j = (int)(*targets);
    if (j<0 || j>=d->int_target) abort();
    a = latent_values[0];
    for (i = 1; i<d->int_target; i++)
    { a = addlogs(a,latent_values[i]);
    }
    l = latent_values[j] - a;
  }

  else if (m->type=='B')
  { l = 0;
    for (j = 0; j<d->N_targets; j++)
    { l += targets[j]!=0 ? -log(1+exp(-latent_values[j])) 
                         : -log(1+exp(latent_values[j]));
    }
  }

  else if (m->type=='R')
  { l = 0;
    for (j = 0; j<d->N_targets; j++)
    { n = noise_variances ? noise_variances[j] : exp (2 * *h->noise[j]);
      a = targets[j] - latent_values[j];
      l -= 0.5*log(2*M_PI*n) + a*a/(2*n);
    }
  }

  else
  { abort();
  }
  
  return l;
}
