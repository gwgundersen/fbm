/* GP-UTIL.C - Utility routines for Gaussian process programs. */

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
#include "gp.h"


/* SET UP REQUIRED RECORD SIZES.  Doesn't set the sizes of the variable-sized
   records. */

void gp_record_sizes
( log_gobbled *logg     /* Structure to hold gobbled data */
)
{
  logg->req_size['P'] = sizeof (gp_spec);
  logg->req_size['M'] = sizeof (model_specification);
  logg->req_size['V'] = sizeof (model_survival);
}


/* REPORT ERROR IF GP SPECS, DATA MODEL, ETC. ARE MISSING. */

void gp_check_specs_present
( gp_spec *gp,			/* GP specifications, or null */
  int need_model,		/* Must model be present? */
  model_specification *m,	/* Model, or null */
  model_survival *v		/* Survival model parameters, or null */
)
{
  if (gp==0)
  { fprintf(stderr,"No specification for Gaussian process model in log file\n");
    exit(1);
  }

  if (m==0 && need_model)
  { fprintf(stderr,"No model specification in log file\n");
    exit(1);
  }

  if (m!=0 && m->type=='V' && v==0)
  { fprintf(stderr,"No hazard specification for survival model\n");
    exit(1);
  }
}


/* RETURN NUMBER OF HYPERPARAMETERS.  Returns the number of hyperparameter
   values for Gaussian process model with the given specification, for
   use with the given data model. */
   
int gp_hyper_count
( gp_spec *gp,			/* Gaussian process specification */
  model_specification *m 	/* Data model */
)
{ 
  int count;
  int l, i;

  count = 0;

  if (gp->has_constant)
  { if (gp->constant.alpha[0]!=0) count += 1;
  }

  if (gp->has_linear) 
  { if (gp->linear.alpha[0]!=0) count += 1;
    if (gp->linear.alpha[1]!=0) 
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->linear_flags[i]&Flag_omit))
        { count += 1;
        }
      }
    }
  }

  if (gp->has_jitter)
  { if (gp->jitter.alpha[0]!=0) count += 1;
  }
  
  for (l = 0; l<gp->N_exp_parts; l++)
  { if (gp->exp[l].scale.alpha[0]!=0) count += 1;
    if (gp->exp[l].relevance.alpha[0]!=0) count += 1;
    if (gp->exp[l].relevance.alpha[1]!=0) 
    { for (i = 0; i<gp->N_inputs; i++)
      { if (!(gp->exp[l].flags[i]&Flag_omit))
        { count += 1;
        }
      }
    }
  }

  if (m!=0 && m->type=='R') 
  { if (m->noise.alpha[0]!=0) count += 1;
    if (m->noise.alpha[1]!=0) count += gp->N_outputs;
  }
  
  return count;
}


/* SET UP POINTERS TO HYPERPARAMETERS.  Sets the pointers in the gp_hypers
   structure to point to the appropriate places in a block of hyperparameters.
   Pointers associated with parts that don't exist are set to zero.  Sometimes
   several pointers point to the same place, when hyperparameters are tied
   together, and sometimes the pointers point to constants, when the prior
   is degenerate.

   The size of the storage block, as returned by the gp_hyper_count 
   procedure, must be stored by the caller in the total_hypers field, and
   the caller must also set the hyper_block field to point to a block of 
   this many double values. */

void gp_hyper_pointers
( gp_hypers *h,			/* Structure to set up pointers in */
  gp_spec *gp,			/* Gaussian process specification */
  model_specification *m 	/* Data model */
)
{ 
  double *b;
  int l, i;

  b = h->hyper_block;
  
  if (!gp->has_constant)
  { h->constant = 0;
  }
  else 
  { if (gp->constant.alpha[0]!=0)
    { h->constant = b++;
    }
    else
    { h->constant = &h->const_constant;
      h->const_constant = log(gp->constant.width);
    }
  }

  if (!gp->has_linear)
  { h->linear_cm = 0;
    for (i = 0; i<gp->N_inputs; i++) h->linear[i] = 0;
  }
  else 
  { if (gp->linear.alpha[0]!=0)
    { h->linear_cm = b++;
    }
    else
    { h->linear_cm = &h->const_linear;
      h->const_linear = log(prior_width_scaled(&gp->linear,gp->N_inputs));
    }
    for (i = 0; i<gp->N_inputs; i++)
    { if (gp->linear_flags[i]&Flag_omit)
      { h->linear[i] = 0;
      }
      else if (gp->linear.alpha[1]!=0)
      { h->linear[i] = b++;
      }
      else
      { h->linear[i] = h->linear_cm;
      }
    }
  }
  
  if (!gp->has_jitter)
  { h->jitter = 0;
  }
  else 
  { if (gp->jitter.alpha[0]!=0)
    { h->jitter = b++;
    }
    else
    { h->jitter = &h->const_jitter;
      h->const_jitter = log(gp->jitter.width);
    }
  }

  for (l = 0; l<gp->N_exp_parts; l++)
  { 
    if (gp->exp[l].scale.alpha[0]!=0)
    { h->exp[l].scale = b++;
    }
    else
    { h->exp[l].scale = &h->exp[l].const_scale;
      h->exp[l].const_scale = log(gp->exp[l].scale.width);
    }

    if (gp->exp[l].relevance.alpha[0]!=0)
    { h->exp[l].rel_cm = b++;
    }
    else
    { h->exp[l].rel_cm = &h->exp[l].const_rel;
      h->exp[l].const_rel = 
        log (prior_width_scaled (&gp->exp[l].relevance, gp->N_inputs));
    }
    for (i = 0; i<gp->N_inputs; i++)
    { if (gp->exp[l].flags[i]&Flag_omit) 
      { h->exp[l].rel[i] = 0;
      }
      else if (gp->exp[l].relevance.alpha[1]==0)
      { h->exp[l].rel[i] = h->exp[l].rel_cm;
      }
      else
      { h->exp[l].rel[i] = b++;
      }
    }
  }

  if (m==0 || m->type!='R')
  { h->noise_cm = 0;
    for (i = 0; i<gp->N_outputs; i++) h->noise[i] = 0;
  }
  else
  { if (m->noise.alpha[0]!=0)
    { h->noise_cm = b++;
    }
    else
    { h->noise_cm = &h->const_noise;
      h->const_noise = log(m->noise.width);
    }
    for (i = 0; i<gp->N_outputs; i++)
    { if (m->noise.alpha[1]!=0)
      { h->noise[i] = b++;
      }
      else
      { h->noise[i] = h->noise_cm;
      }
    }
  }
}
