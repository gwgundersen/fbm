/* MODEL.C - Procedures dealing with modeling target value. */

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
#include "data.h"
#include "model.h"
#include "rand.h"


/* CHECK THAT OUTPUT VALUES ARE APPROPRIATE FOR THE MODEL, ETC.  
   An error message is produced and the program halted if not.  No model 
   at all, indicated by a null pointer, is treated the same an 'R' model. 
   The last argument is a string with the model IDs for which the log
   of sum of exponentials scheme is allowed, or null if it's always allowed. 

   The number of values needed after log sum exp computations is returned
   as the result (but this may be ignored, of course). */

int model_values_check
( model_specification *m,	/* Model specification */
  data_specifications *d,	/* Data source specification */
  int n_values,			/* Number of output values in model */
  char *log_sum_exp_ok		/* Models for which log sum exp scheme is ok */
)
{ 
  int t, n;

  if (m!=0 && m->type=='B' && d->int_target!=2)
  { fprintf(stderr,"Data for binary targets must be specified to be binary\n");
    exit(1);
  }

  if (m!=0 && m->type=='C' && d->int_target==0)
  { fprintf(stderr,"Data for class must have specified number of values\n");
    exit(1);
  }

  if (m!=0 && m->type=='N' && d->int_target==0)
  { fprintf(stderr,"Data for count targets must be specified to be integer\n");
    exit(1);
  }

  t = m==0 || m->type==0 ? 'R' : m->type;

  switch (t)
  { 
    case 'B': case 'N': case 'R': case 'r':
    { n = d->N_targets;
      break;
    }

    case 'C': 
    { n = d->int_target;
      break;
    }

    case 'V': 
    { n = 1;
      break;
    }

    default:
    { fprintf(stderr,"Unknown data model: %c\n",t);
      exit(1);
    }
  }

  if (n_values!=n)
  { if (log_sum_exp_ok!=0 && strchr(log_sum_exp_ok,t)==0 || n_values%n!=0)
    { fprintf(stderr,
       "Number of model outputs is not appropriate for data model (%d %d %c %s)\n",
       n_values, n, t, log_sum_exp_ok);
      exit(1);
    }
  }

  return n;  
}


/* RANDOMLY GENERATE TARGETS ACCORDING TO THE MODEL.  No model is treated
   like an 'R' model with zero noise. */

void model_gen
( model_specification *m,	/* Model specification */
  data_specifications *d,	/* Data source specification */
  double *noise,		/* Noise variances, 0 if targets are not real */
  double *values,		/* Values used to generate targets */
  int n_values,			/* Number of such values */
  double *targets		/* Place to store generated targets */
)
{
  double *t;
  double mx;
  int j, k;

  t = model_values(m,d,values,n_values);

  switch (m==0 ? 0 : m->type)
  {
    case 0:
    { for (j = 0; j<d->N_targets; j++)
      { targets[j] = t[j];
      }
      break;
    }
    case 'R': 
    { for (j = 0; j<d->N_targets; j++)
      { targets[j] = t[j] + rand_gaussian() * sqrt(noise[j]);
      }
      break;
    }
    case 'B':
    { for (j = 0; j<d->N_targets; j++)
      { targets[j] = rand_uniform() < 1/(1+exp(-t[j]));
      }
      break;
    }
    case 'r': 
    { for (j = 0; j<d->N_targets-1; j++)
      { targets[j] = t[j] + rand_gaussian() * sqrt(noise[j]);
      }
      targets[j] = rand_uniform() < 1/(1+exp(-t[j]));
      break;
    }
    case 'N':
    { for (j = 0; j<d->N_targets; j++)
      { targets[j] = rand_poisson(exp(t[j]));
      }
      break;
    }
    case 'C':
    { 
      if (d->N_targets!=1) abort();

      mx = t[0];
      for (j = 1; j<d->int_target; j++)
      { if (t[j]>mx)   
        { mx = t[j];
        }
      }

      for (j = 0; j<d->int_target; j++)
      { t[j] = exp(t[j]-mx);
      }

      *targets = rand_pickd(t,d->int_target);

      break;
    }

    case 'V':
    { fprintf(stderr,"Can't generate targets for survival models\n");
      exit(1);
    }

    default:
    { fprintf(stderr,"Unknown data model: %c\n",m->type);
      exit(1);
    }
  }
}


/* COMPUTE LOG LIKELIHOOD FOR TARGETS GIVEN LATENT VALUES FOR ONE CASE. */

double model_likelihood
( model_specification *m,	/* Specification for model used */
  data_specifications *d,	/* Specification for data */
  double *noise,		/* Noise variances, 0 if targets are not real */
  double *values,		/* Values used to generate targets */
  int n_values,			/* Number of such values */
  double *targets		/* Target values being modeled */
)
{ 
  double *t;
  double l, a, s;
  int i, j, k;

  t = model_values(m,d,values,n_values);

  if (m->type=='C')
  { j = (int)(*targets);
    if (j<0 || j>=d->int_target) abort();
    a = t[0];
    for (i = 1; i<d->int_target; i++)
    { a = addlogs(a,t[i]);
    }
    l = t[j] - a;
  }

  else if (m->type=='B')
  { l = 0;
    for (j = 0; j<d->N_targets; j++)
    { l += targets[j]!=0 ? -log(1+exp(-t[j])) 
                         : -log(1+exp(t[j]));
    }
  }

  else if (m->type=='N')
  { l = 0;
    for (j = 0; j<d->N_targets; j++)
    { if (targets[j]>=0)
      { l += targets[j]*t[j] - lgamma(targets[j]+1) 
              - exp(t[j]);
      }
      else
      { s = - exp(t[j]);
        for (i = 1; i <= -targets[j]; i++)
        { s = addlogs (s, i*t[j]-lgamma(i+1)-exp(t[j]));
        }
        l += s;
      }
    }
  }

  else if (m->type=='R')
  { l = 0;
    for (j = 0; j<d->N_targets; j++)
    { a = targets[j] - t[j];
      l -= 0.5*log(2*M_PI*noise[j]) + a*a/(2*noise[j]);
    }
  }

  else
  { abort();
  }

  return l;
}


/* COMPUTE VALUES USED BY MODEL FROM OUTPUTS.  Returns a pointer to these
   values, which may be just the output vector, or may be a statically
   allocated area which will be overwritten by the next call of this
   procedure. */

double *model_values
( model_specification *m,	/* Specification for model used */
  data_specifications *d,	/* Specification for data */
  double *values,		/* Values used to generate targets */
  int n_values			/* Number of such values */
)
{
  static double *t = 0;
  static int tn;
  int j, k, n, q;
  double mx;

  /* See how many values we need for finding likelihood of targets. */

  n = model_values_check(m,d,n_values,0);

  /* Just return outputs if they match required values. */

  if (n==n_values) 
  { return values;
  }

  /* Allocate space the first time. */

  if (t==0 || tn!=n)
  { if (t!=0) free (t);
    t = chk_alloc(n,sizeof(double));
    tn = n;
  }

  /* Obtain these values using the log-sum-exp scheme. */

  q = n_values/n;

  for (j = 0; j<n; j++)
  { 
    mx = values[j*q];
    for (k = 1; k<q; k++)
    { if (values[j*q+k]>mx)
      { mx = values[j*q+k];
      }
    }

    t[j] = 0;
    for (k = 0; k<q; k++)
    { t[j] += exp(values[j*q+k]-mx);
    }

    t[j] = log(t[j]) + mx;
  }

  return t;
}
