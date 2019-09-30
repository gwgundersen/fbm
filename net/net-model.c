/* NET-MODEL.C - Module dealing with the interpretation of network outputs. */

/* Copyright (c) 1995 by Radford M. Neal 
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
#include "net.h"
#include "rand.h"


/* RETURN NUMBER OF TARGET VALUES REQUIRED.  No model at all is treated
   the same a the 'R' model. */

int net_model_targets
( net_arch *a		/* Network architecture */
)
{ switch (a->data_model)
  { 
    case 'B': case 'R': case 0: return a->N_outputs;

    case 'C': return 1;

    default:
    { fprintf(stderr,"Unknown data model: %c\n",a->data_model);
      exit(1);
    }
  }
}


/* COMPUTE LOG PROBABILITY OF TARGETS AND/OR ITS DERIVATIVES.  Computes the 
   log of the probability (or probability density) of the observed target 
   values given the observed inputs, as defined by the current network outputs 
   for those inputs, and the derivative of minus the log probability with 
   respect to the network outputs. 

   If the 'op' parameter is 2, the probability may be computed ignoring factors 
   that depend only on the overall and per-unit noise hyperparameters, and this
   freedom will be used to ensure that values greater than zero will not be
   returned; if 'op' is 1, factors that are constant or depend only on the 
   alphas may be ignored; if 'op' is 0, the exact probability will be computed.

   By passing zero pointers, computation of either the log probability or of
   its derivatives may be suppressed, with a possible saving in time. 
 
   A data model must be specified to use this procedure. */

void net_model_prob
( net_values *v,	/* Values for units in network */
  double *t,		/* Target values */
  double *pr,		/* Place to store log probability, zero if not wanted */
  net_values *dp,	/* Place to store log probability derivatives, or zero*/
  net_arch *a,		/* Network architecture */
  net_priors *p,	/* Network priors, including noise priors */
  net_sigmas *s,	/* Hyperparameters, including noise sigmas */
  int op		/* Can we ignore some factors? */
)
{
  extern double lgamma(double);

  double p1, d, x, alpha, z;
  int i;

  switch (a->data_model)
  {
    case 'B':  /* Binary data values */
    { 
      if (pr) *pr = 0;

      for (i = 0; i<a->N_outputs; i++)
      { p1 = 1 / (1+exp(-v->o[i]));
        if (pr) 
        { *pr += t[i]==0 ? log(1-p1) : log(p1);
        }
        if (dp)
        { dp->o[i] = p1 - t[i];
        }
      }

      break;
    }

    case 'C':  /* Single class with multiple possible values */
    {
      z = v->o[0];

      for (i = 1; i<a->N_outputs; i++)
      { z = addlogs(z,v->o[i]);
      }

      if (pr) 
      { *pr = v->o[(int)*t] - z;
      }

      if (dp)
      { for (i = 0; i<a->N_outputs; i++)
        { dp->o[i] = exp (v->o[i] - z);
        }
        dp->o[(int)*t] -= 1;
      }

      break;
    }
  
    case 'R':  /* Real-valued target */
    { 
      alpha = p->noise.alpha[2];

      if (alpha==0) /* Gaussian distribution for noise */
      { 
        if (pr) *pr = 0;

        for (i = 0; i<a->N_outputs; i++)
        { d = (v->o[i] - t[i]) / s->noise[i];
          if (d<-1e10) d = -1e10;
          if (d>+1e10) d = +1e10;
          if (pr) 
          { *pr -= 0.5*d*d;
            if (op<2) *pr -= log(s->noise[i]);
          }
          if (dp)
          { dp->o[i] = d / s->noise[i];
          }
        }

        if (pr && op<1) 
        { *pr -= 0.5 * a->N_outputs * log(2*M_PI);
        }
      }

      else /* Student t distribution for noise */
      {
        if (pr) *pr = 0;

        for (i = 0; i<a->N_outputs; i++)
        { d = (v->o[i] - t[i]) / s->noise[i];
          if (d<-1e10) d = -1e10;
          if (d>+1e10) d = +1e10;
          x = 1 + d*d/alpha;
          if (pr) 
          { *pr -= ((alpha+1)/2) * log(x);
            if (op<2) 
            { *pr -= log(s->noise[i]);
            }
          }
          if (dp)
          { dp->o[i] = ((alpha+1)/alpha) * (d/s->noise[i]) / x;
          }
        }

        if (pr && op<1) 
        { *pr += a->N_outputs * (lgamma((alpha+1)/2) - lgamma(alpha/2) 
                                  - 0.5*log(M_PI*alpha));
        }
      }

      break;
    }
    
    case 0:
    { fprintf(stderr,"Network has no data model defined\n");
      exit(1);
    }

    default:
    { fprintf(stderr,"Unknown data model: %c\n",a->data_model);
      exit(1);
    }
  }
}


/* COMPUTE MAXIMUM LOG LIKELIHOOD SECOND DERIVATIVES.  Computes the maximum
   values of the second derivatives of minus the log probability of the targets
   in the training set with respect to the outputs of the net, for the current 
   values of the hyperparameters.  The maximum is with respect to possible
   targets and possible values of the outputs of the net.
 
   A data model must be specified to use this procedure. */

void net_model_max_second
( net_value *msd,	/* Place to store maximum second derivatives */
  net_arch *a,		/* Network architecture */
  net_priors *p,	/* Network priors, including noise priors */
  net_sigmas *s		/* Hyperparameters, including noise sigmas */
)
{
  double alpha;
  int i;

  switch (a->data_model)
  {
    case 'B':  /* Binary data values */
    {
      for (i = 0; i<a->N_outputs; i++)
      { msd[i] = 0.25;
      }

      break;
    }
    
    case 'C':  /* Single class with multiple possible values */
    {
      for (i = 0; i<a->N_outputs; i++)
      { msd[i] = 0.25;
      }

      break;
    }

    case 'R':  /* Real-valued target */
    {
      alpha = p->noise.alpha[2];

      for (i = 0; i<a->N_outputs; i++)
      { msd[i] = alpha==0 ? 1 / (s->noise[i] * s->noise[i])
                          : (alpha+1) / (alpha * s->noise[i] * s->noise[i]);
      }

      break;
    }

    case 0:
    { fprintf(stderr,"Network has no data model defined\n");
      exit(1);
    }

    default:
    { fprintf(stderr,"Unknown data model: %c\n",a->data_model);
      exit(1);
    }
  }
}


/* MAKE GUESSES AT TARGET VALUES.  The guesses are based on the outputs
   of the network (which must already have been computed), and possibly
   on the noise sigma/prior.  Guesses can be either "mean" values, or
   values randomly drawn from the target distribution, depending on the
   setting of the 'rnd' argument.

   For binary targets, the "mean" guess is the probability of the target 
   being one.  The random guess is a 0 or 1, with 1 having the specified
   probability.

   For class targets, the "mean" guess is a vector of probabilities for the 
   various classes.  (Note that this requires an array for the target, not 
   a single value.)  The random guess is one of the classes (numbered from
   zero on up), chosen according to this distribution.  (In this case, only
   a single target is produced.)
 
   For real-valued data, the "mean" guess is just the network output.  The
   random guess is chosen from a Gaussian distribution with this mean and
   with a standard deviation determined from the noise hyperparameters.
   (When each case has a different noise level, the noise is effectively 
   from a t distribution.) 

   If no data model is specified, a real-valued data model with no noise
   is assumed. */

double net_model_guess
( net_values *v,	/* Values for units in network */
  double *t,		/* Place to store guesses at targets */
  net_arch *a,		/* Network architecture */
  net_priors *p,	/* Network priors, including noise priors */
  net_sigmas *s,	/* Hyperparameters, including noise sigmas */
  int rnd		/* Make random guess, rather than use mean? */
)
{
  double z, pr, r, noise, alpha;
  int i;

  switch (a->data_model)
  {
    case 'B':  /* Binary data values */
    { 
      for (i = 0; i<a->N_outputs; i++)
      { pr = 1 / (1+exp(-v->o[i]));
        t[i] = rnd ? rand_uniform()<pr : pr;
      }

      break;
    }

    case 'C':  /* Single class with multiple possible values */
    {
      z = v->o[0];

      for (i = 1; i<a->N_outputs; i++)
      { z = addlogs(z,v->o[i]);
      }

      if (rnd)
      { r = rand_uniform();
        *t = 0; /* just in case */
      }

      for (i = 0; i<a->N_outputs; i++)
      { pr = exp (v->o[i] - z);
        if (rnd)
        { r -= pr;
          if (r<0) 
          { *t = i;
            break;
          }
        }
        else
        { t[i] = pr;
        }
      }

      break;
    }
  
    case 'R': case 0:  /* Real-valued target */
    { 
      for (i = 0; i<a->N_outputs; i++)
      { t[i] = v->o[i];
        if (rnd && a->data_model!=0)
        { noise = s->noise[i];
          alpha = p->noise.alpha[2];
          if (alpha!=0)
          { noise /= sqrt (rand_gamma(alpha/2) / (alpha/2));
          }
          t[i] += noise * rand_gaussian();
        }
      }

      break;
    }

    default:
    { fprintf(stderr,"Unknown data model: %c\n",a->data_model);
      exit(1);
    }
  }
}
