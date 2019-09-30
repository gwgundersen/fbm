/* GP-PRIOR.C - Routines dealing with priors for Gaussian processes. */

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


/* GENERATE VALUES FOR HYPERPARAMETERS. */

void gp_prior_generate(
  gp_hypers *h,		/* Place to store generated hyperparameters */
  gp_spec *gp,		/* Specification for Gaussian process */
  model_specification *m, /* Specification for data model */
  int fix,		/* Should parameters be fixed, rather than random? */
  double scale_value,	/* Value to fix scale parameters at, 0 if unspecified */
  double relevance_value /* Value to fix relevance parameters at */
)
{
  int i, l;

  /* Pick Gaussian process hyperparameters. */

  if (gp->has_constant)
  { if (gp->constant.alpha[0]!=0)
    { *h->constant = log (
        !fix ? prior_pick_sigma(gp->constant.width,gp->constant.alpha[0])
             : scale_value==0 ? gp->constant.width : scale_value);
    }
  }

  if (gp->has_linear)
  { 
    if (gp->linear.alpha[0]!=0)
    { *h->linear_cm = log (
        !fix ? prior_pick_sigma(prior_width_scaled(&gp->linear,gp->N_inputs),
                                gp->linear.alpha[0])
             : scale_value==0 ? gp->linear.width : scale_value*relevance_value);
    }

    if (gp->linear.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->linear_flags[i]&Flag_omit))
        { *h->linear[i] = 
           log (!fix ? prior_pick_sigma(exp(*h->linear_cm),gp->linear.alpha[1])
           : scale_value==0 ? exp(*h->linear_cm) : scale_value*relevance_value);
        }
      }
    }
  }

  if (gp->has_jitter)
  { if (gp->jitter.alpha[0]!=0)
    { *h->jitter = log (
        !fix ? prior_pick_sigma(gp->jitter.width,gp->jitter.alpha[0])
             : gp->jitter.width);
    }
  }
 
  for (l = 0; l<gp->N_exp_parts; l++)
  {   
    if (gp->exp[l].scale.alpha[0]!=0)
    { *h->exp[l].scale = log (
       !fix ? prior_pick_sigma(gp->exp[l].scale.width,gp->exp[l].scale.alpha[0])
            : scale_value==0 ? gp->exp[l].scale.width : scale_value);
    }

    if (gp->exp[l].relevance.alpha[0]!=0)
    { *h->exp[l].rel_cm = log (
       !fix ? prior_pick_sigma(prior_width_scaled(&gp->exp[l].relevance,
                                                  gp->N_inputs),
                               gp->exp[l].relevance.alpha[0])
            : scale_value==0 ? prior_width_scaled(&gp->exp[l].relevance,
                                                  gp->N_inputs)
                             : relevance_value);
    }

    if (gp->exp[l].relevance.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->exp[l].flags[i]&Flag_omit))
        { *h->exp[l].rel[i] = log (
           !fix ? prior_pick_sigma (exp(*h->exp[l].rel_cm),
                                    gp->exp[l].relevance.alpha[1])
                : scale_value==0 ? exp(*h->exp[l].rel_cm) : relevance_value);
        }
      }
    }
  }

  /* Pick hyperparameters pertaining to data model. */

  if (m!=0 && m->type=='R') 
  {
    if (m->noise.alpha[0]!=0)
    { *h->noise_cm = log ( 
       !fix ? prior_pick_sigma(m->noise.width,m->noise.alpha[0])
            : scale_value==0 ? m->noise.width : scale_value);
    }

    if (m->noise.alpha[1]!=0)
    { for (i = 0; i<gp->N_outputs; i++)
      { *h->noise[i] = log (
         !fix ? prior_pick_sigma(exp(*h->noise_cm),m->noise.alpha[1])
              : scale_value==0 ? exp(*h->noise_cm) : scale_value);
      }
    }
  }
}


/* FIND LOG PRIOR DENSITY FOR HYPERPARAMETERS.  The log density found is for 
   the hyperparameters in logarithmic form (ie, as logs of the 'width' form).
   The density omits the degenerate contributions of hyperparameters tied
   to being exactly equal to higher-level hyperparameters (or to constants). 

   The density computed may omit constant factors (depending only one the
   shape parameters) if the last parameter is zero; the exact result is 
   returned if the last parameter is non-zero. */

double gp_log_prior
( gp_hypers *h,			/* Hyperparameter values */
  gp_spec *gp,			/* Specification for Gaussian process */
  model_specification *m,	/* Specification for data model (or zero) */
  int ex			/* Is exact value required? */
)
{ 
  double lp;
  int l, i;

  lp = 0;

  if (gp->has_constant)
  { if (gp->constant.alpha[0]!=0) 
    { lp += gp_gdens (gp->constant.alpha[0], gp->constant.width, 
                      *h->constant, ex);
    }
  }

  if (gp->has_linear)
  { if (gp->linear.alpha[0]!=0)
    { lp += gp_gdens (gp->linear.alpha[0],
              prior_width_scaled(&gp->linear,gp->N_inputs), *h->linear_cm, ex);
    }
    if (gp->linear.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->linear_flags[i]&Flag_omit))
        { lp += gp_gdens (gp->linear.alpha[1], exp(*h->linear_cm), 
                          *h->linear[i], ex);
        }
      }
    }
  }

  if (gp->has_jitter)
  { if (gp->jitter.alpha[0]!=0) 
    { lp += gp_gdens (gp->jitter.alpha[0], gp->jitter.width, *h->jitter, ex);
    }
  }

  for (l = 0; l<gp->N_exp_parts; l++)
  { 
    if (gp->exp[l].scale.alpha[0]!=0) 
    { lp += gp_gdens (gp->exp[l].scale.alpha[0], gp->exp[l].scale.width, 
                      *h->exp[l].scale, ex);
    }

    if (gp->exp[l].relevance.alpha[0]!=0) 
    { lp += gp_gdens (gp->exp[l].relevance.alpha[0], 
                      prior_width_scaled(&gp->exp[l].relevance,gp->N_inputs),
                      *h->exp[l].rel_cm, ex);
    }

    if (gp->exp[l].relevance.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->exp[l].flags[i]&Flag_omit))
        { lp += gp_gdens (gp->exp[l].relevance.alpha[1], exp(*h->exp[l].rel_cm),
                          *h->exp[l].rel[i], ex);
        }
      }
    }
  }

  if (m!=0 && m->type=='R')
  { 
    if (m->noise.alpha[0]!=0)
    { lp += gp_gdens (m->noise.alpha[0], m->noise.width, *h->noise_cm, ex);
    }

    if (m->noise.alpha[1]!=0)
    { for (i = 0; i<gp->N_outputs; i++)
      { lp += gp_gdens (m->noise.alpha[1], exp(*h->noise_cm), *h->noise[i], ex);
      }
    }
  }

  return lp;
}


/* FIND GRADIENT OF MINUS LOG PRIOR DENSITY.  The derivatives of minus the
   log density of the hyperparameters (in logarithmetic form) are found,
   and stored in the last parameter. 

   The derivatives are found by first setting them all to zero, and then
   looking at each contribution to the prior density (corresponding to one
   conditional probability of gamma form), and adding to the derivatives
   of the (up to) two hyperparameters that it may involve. */

void gp_prior_grad
( gp_hypers *h,			/* Hyperparameter values */
  gp_spec *gp,			/* Specification for Gaussian process */
  model_specification *m,	/* Specification for data model (or zero) */
  gp_hypers *d			/* Place to store derivatives */
)
{ 
  int l, i;

  for (i = 0; i<d->total_hypers; i++) d->hyper_block[i] = 0;

  if (gp->has_constant)
  { if (gp->constant.alpha[0]!=0) 
    { gp_gdiff (gp->constant.alpha[0], gp->constant.width, *h->constant, 
                0, d->constant);
    }
  }

  if (gp->has_linear)
  { if (gp->linear.alpha[0]!=0)
    { gp_gdiff (gp->linear.alpha[0],
                prior_width_scaled(&gp->linear,gp->N_inputs), *h->linear_cm, 
                0, d->linear_cm);
    }
    if (gp->linear.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->linear_flags[i]&Flag_omit))
        { gp_gdiff (gp->linear.alpha[1], exp(*h->linear_cm), *h->linear[i], 
                    gp->linear.alpha[0]==0 ? 0 : d->linear_cm, d->linear[i]);
        }
      }
    }
  }

  if (gp->has_jitter)
  { if (gp->jitter.alpha[0]!=0) 
    { gp_gdiff (gp->jitter.alpha[0], gp->jitter.width, *h->jitter, 
                0, d->jitter);
    }
  }

  for (l = 0; l<gp->N_exp_parts; l++)
  { 
    if (gp->exp[l].scale.alpha[0]!=0) 
    { gp_gdiff (gp->exp[l].scale.alpha[0], gp->exp[l].scale.width, 
                *h->exp[l].scale, 0, d->exp[l].scale);
    }

    if (gp->exp[l].relevance.alpha[0]!=0) 
    { gp_gdiff (gp->exp[l].relevance.alpha[0], 
                prior_width_scaled(&gp->exp[l].relevance,gp->N_inputs),
                *h->exp[l].rel_cm,
                0, d->exp[l].rel_cm);
    }

    if (gp->exp[l].relevance.alpha[1]!=0)
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->exp[l].flags[i]&Flag_omit))
        { gp_gdiff (gp->exp[l].relevance.alpha[1], 
                    exp(*h->exp[l].rel_cm), 
                    *h->exp[l].rel[i], 
                    gp->exp[l].relevance.alpha[0]==0 ? 0 : d->exp[l].rel_cm,
                    d->exp[l].rel[i]);
        }
      }
    }
  }

  if (m!=0 && m->type=='R')
  { 
    if (m->noise.alpha[0]!=0)
    { gp_gdiff (m->noise.alpha[0], m->noise.width, *h->noise_cm, 
                0, d->noise_cm);
    }

    if (m->noise.alpha[1]!=0)
    { for (i = 0; i<gp->N_outputs; i++)
      { gp_gdiff (m->noise.alpha[1], exp(*h->noise_cm), *h->noise[i],
                  m->noise.alpha[0]==0 ? 0 : d->noise_cm, d->noise[i]);
      }
    }
  }
}


/* FIND LOG DENSITY FOR HYPERPARAMETER.  Computes the log of the density for 
   the log of a hyperparameter which in precision form has a gamma density with
   shape and scale parameter as given.  Constant factors (depending only on
   the shape parameter) may be omitted if the last parameter is zero. 

   The log density for a gamma variable, x, is given by

       (alpha/2 - 1) * log(x) - x * alpha/(2*mean) 
          + (alpha/2) * log(alpha/(2*mean)) - lgamma(alpha/2) 

   where lgamma is the log of the gamma function.  The last terms depends
   only on alpha, and so can be omitted if the last parameter is zero.

   It is in their "precision" form that the hyperparameters have gamma 
   prior densities, as above.  The value actually represented is the log 
   of the "width" form, that is v = log(1/sqrt(x)) = -(1/2)*log(x).  The
   log density for v is therefore the log density for x shown above plus
   log(x) + log(2).  The mean precision in terms of the 'width' in the 
   prior specification is 1/(width*width).
*/

double gp_gdens
( double alpha,
  double width,
  double v,
  int exact
)
{ 
  double mean, lp;
  extern double lgamma(double); 

  mean = 1/(width*width);
 
  lp = (alpha/2) * (-2*v) - exp(-2*v) * alpha/(2*mean) 
         + (alpha/2) * log(alpha/(2*mean));

  if (exact)
  { lp += log(2.0) - lgamma(alpha/2);
  }

  return lp;
}


/* FIND DERIVATIVES OF MINUS LOG DENSITY.  Finds the derivatives of minus
   the log prior density for a hyperparameter in logarithmic form with respect 
   to the log of the width parameter and the hyperparameter itself (in log
   form).  These derivatives are added to the values to which pointers are
   passed.  A zero pointer suppresses calculation of the corresponding 
   derivative. */

void gp_gdiff
( double alpha,
  double width,
  double v,
  double *dwidth,
  double *dv
)
{
  double mean;

  mean = 1/(width*width);

  if (dwidth!=0)
  { *dwidth += exp(-2*v)*alpha/mean - alpha;
  }

  if (dv!=0)
  { *dv += alpha - exp(-2*v)*alpha/mean;
  }
}
