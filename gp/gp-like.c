/* GP-LIKE.C - Likelihood computations for Gaussian Process models. */  
 
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
  double l, n, a, s;
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

  else if (m->type=='N')
  { l = 0;
    for (j = 0; j<d->N_targets; j++)
    { if (targets[j]>=0)
      { l += targets[j]*latent_values[j] - lgamma(targets[j]+1) 
              - exp(latent_values[j]);
      }
      else
      { s = - exp(latent_values[j]);
        for (i = 1; i <= -targets[j]; i++)
        { s = addlogs (s, i*latent_values[j]-lgamma(i+1)-exp(latent_values[j]));
        }
        l += s;
      }
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
