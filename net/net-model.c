/* NET-MODEL.C - Module dealing with the interpretation of network outputs. */

/* Copyright (c) 1995, 1996, 1999 by Radford M. Neal 
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
#include "net.h"
#include "rand.h"


/* CONSTANTS INVOLVING PI. */

#ifndef M_PI
#define M_PI 3.14159265358979323846	/* Define pi, if not defined already */
#endif

#define Log2pi  1.83787706640934548356	/* Log(2*M_PI) */


/* COMPUTE LOG PROBABILITY OF TARGETS AND/OR ITS DERIVATIVES.  Computes the 
   log of the probability (or probability density) of the observed target 
   values given the observed inputs, as defined by the current network outputs 
   for those inputs, and the derivative of minus the log probability with 
   respect to the network outputs. 
 
   A data model must be specified to use this procedure.  For survival models
   with piecewise constant hazard, this procedure should be called several
   times, once for each piece, with a target value that is the survival
   time after the start of the piece (negative if censored).  The log
   probabilities from these calls should then be added together (though the 
   derivatives will generally be required in separate form).

   If the 'op' parameter is 2, the probability may be computed ignoring factors 
   that depend only on the overall and per-unit noise hyperparameters, and this
   freedom will be used to ensure that values greater than zero will not be
   returned, except for survival models; if 'op' is 1, factors that are 
   constant or depend only on the alphas may be ignored; if 'op' is 0, the 
   exact probability will be computed.

   By passing zero pointers, computation of either the log probability or of
   its derivatives may be suppressed, with a possible saving in time. */

void net_model_prob
( net_values *v,	/* Values for units in network */
  double *t,		/* Target values, fudged for piecewise const hazard */
  double *pr,		/* Place to store log probability, zero if not wanted */
  net_values *dp,	/* Place to store log probability derivatives, or zero*/
  net_arch *a,		/* Network architecture */
  model_specification *m, /* Data model */
  model_survival *sv,	/* Type of hazard function for survival model, or null*/
  net_sigmas *s,	/* Hyperparameters, including noise sigmas */
  int op		/* Can we ignore some factors? */
)
{
  extern double lgamma(double);

  static double alpha_saved=0.0;/* Constant is already computed for this alpha*/
  static double cnst;		/* Saved value of this constant */

  double d, x, alpha, z;
  int i;

  switch (m->type)
  {
    case 'B':  /* Binary data values */
    { 
      if (pr) *pr = 0;

      for (i = 0; i<a->N_outputs; i++)
      { double oi;
        oi = v->o[i];
        if (t[i]==0)
        { if (oi<0)
          { if (pr) *pr -= log(1+exp(oi));
            if (dp) dp->o[i] = 1 - 1/(1+exp(oi));
          }
          else
          { if (pr) *pr -= oi + log(1+exp(-oi));
            if (dp) dp->o[i] = 1/(1+exp(-oi));
          }
        }
        else
        { if (oi<0)
          { if (pr) *pr -= -oi + log(1+exp(oi));
            if (dp) dp->o[i] = -1/(1+exp(oi));
          }
          else
          { if (pr) *pr -= log(1+exp(-oi));
            if (dp) dp->o[i] = -1 + 1/(1+exp(-oi));
          }
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
      alpha = m->noise.alpha[2];

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
        { *pr -= 0.5 * a->N_outputs * Log2pi;
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
        { if (alpha!=alpha_saved)
          { cnst = lgamma((alpha+1)/2) - lgamma(alpha/2) - 0.5*log(M_PI*alpha);
            alpha_saved = alpha;
          }
          *pr += a->N_outputs * cnst;
        }
      }

      break;
    }

    case 'V':  /* Survival model */
    { 
      int censored;
      double m, ot, ho;

      if (t[0]<0)
      { censored = 1;
        ot = -t[0];
      }
      else
      { censored = 0;
        ot = t[0];
      }

      m = v->o[0]+log(ot);

      if (m>200.0)
      { ho = exp(200.0);
        if (pr) *pr = - ho;
        if (dp) dp->o[0] = 0;
      }
      else
      { 
        ho = exp(m);

        if (pr) *pr = - ho;
        if (dp) dp->o[0] = ho;

        if (!censored)
        { if (pr) *pr += v->o[0];
          if (dp) dp->o[0] -= 1;
        }
      }

      break;
    }
    
    case 0:
    { fprintf(stderr,"Network has no data model defined\n");
      exit(1);
    }

    default:
    { fprintf(stderr,"Unknown data model: %c\n",m->type);
      exit(1);
    }
  }
}


/* COMPUTE MAXIMUM LOG LIKELIHOOD SECOND DERIVATIVES.  Computes the maximum
   values of the second derivatives of minus the log probability of the targets
   in the training set with respect to the outputs of the net, for the current 
   values of the hyperparameters.  The maximum is with respect to possible 
   values of the outputs of the net (the real outputs must not be looked at,
   or validity of the Markov chain methods would be undermined), and of the
   true targets (which could be looked at, but which aren't here).  In some 
   cases, one must make do with approximations.
 
   A data model must be specified to use this procedure. */

void net_model_max_second
( net_value *msd,	/* Place to store maximum second derivatives */
  net_arch *a,		/* Network architecture */
  model_specification *m, /* Data model */
  model_survival *sv,	/* Type of hazard function for survival model, or null*/
  net_sigmas *s		/* Hyperparameters, including noise sigmas */
)
{
  double alpha;
  int i;

  switch (m->type)
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
      alpha = m->noise.alpha[2];

      for (i = 0; i<a->N_outputs; i++)
      { msd[i] = alpha==0 ? 1 / (s->noise[i] * s->noise[i])
                          : (alpha+1) / (alpha * s->noise[i] * s->noise[i]);
      }

      break;
    }

    case 'V':  /* Survival data */
    {
      msd[0] = 1; /* Rather crude, but not completely arbitrary */

      break;
    }

    case 0:
    { fprintf(stderr,"Network has no data model defined\n");
      exit(1);
    }

    default:
    { fprintf(stderr,"Unknown data model: %c\n",m->type);
      exit(1);
    }
  }
}


/* MAKE GUESSES AT TARGET VALUES.  The guesses are based on the outputs
   of the network (which must already have been computed), except for
   piecewise-constant survival models (see below).  The noise sigma/prior
   may also be relevant.  Guesses can be either "mean" values, or values 
   randomly drawn from the target distribution, depending on the setting 
   of the 'rnd' argument.

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

   For a survival model, the "mean" guess is the mean of the distribution
   of survival times, the random guess is a value drawn from this distribution.
   For piecewise-constant hazards, the network output at the start of this
   procedure is ignored.  The output is instead evaluated several times
   by this procedure, with different settings for the first input (time).

   If no data model is specified, a real-valued data model with no noise
   is assumed. */

void net_model_guess
( net_values *v,	/* Values for units in network */
  double *t,		/* Place to store guesses at targets */
  net_arch *a,		/* Network architecture */
  model_specification *m, /* Data model */
  model_survival *sv,	/* Type of hazard function for survival model, or null*/
  net_params *params,	/* Network parameters (used only for pw-const-hazard) */
  net_sigmas *s,	/* Hyperparameters, including noise sigmas */
  int rnd		/* Make random guess, rather than use mean? */
)
{
  double z, pr, r, noise, alpha;
  int i;

  switch (m->type)
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
        if (rnd && m->type!=0)
        { noise = s->noise[i];
          alpha = m->noise.alpha[2];
          if (alpha!=0)
          { noise /= sqrt (rand_gamma(alpha/2) / (alpha/2));
          }
          t[i] += noise * rand_gaussian();
        }
      }

      break;
    }
  
    case 'V': 
    {
      switch (sv->hazard_type)
      {
        case 'C':  /* Constant hazard */
        {
          double h;
          h = exp(v->o[0]);
          t[0] = rnd ? rand_exp()/h : 1/h;
          break;
        }

        case 'P':  /* Piecewise constant hazard */ 
        {
          double t0, t1, h, pr;
          int w;

          t0 = 0;
          t1 = sv->time[0];
          v->i[0] = sv->log_time ? log(t1) : t1;
 
          t[0] = 0;
          pr = 1;
          w = 0;

          for (;;)
          {
            net_func (v, 0, a, params);
            h = exp(v->o[0]);
            
            if (rnd)
            {
              t[0] = t0 + rand_exp()/h;
              if (t1==-1 || t[0]<=t1) break;
            }
            else
            {
              t[0] += pr * (t0 + 1/h);
              if (t1==-1) 
              { break;
              }
              pr *= exp(-h*(t1-t0));
              t[0] -= pr * (t1 + 1/h);
            }

            t0 = t1;
            w += 1;
          
            if (sv->time[w]==0) 
            { t1 = -1;
              v->i[0] = sv->log_time ? log(t0) : t0;
            }
            else
            { t1 = sv->time[w];
              v->i[0] = sv->log_time ? (log(t0)+log(t1))/2 : (t0+t1)/2;
            }
          }

          break;
        }

        default: 
        { fprintf(stderr,"Unknown hazard type: %c\n",sv->hazard_type);
          exit(1);
        }
      }

      break;
    }

    default:
    { fprintf(stderr,"Unknown data model: %c\n",m->type);
      exit(1);
    }
  }
}
