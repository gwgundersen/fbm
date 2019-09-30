/* NET-UTIL.C - Various utility procedures for use in neural network code. */

/* Copyright (c) 1995, 1996, 2001 by Radford M. Neal 
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
#include "data.h"
#include "prior.h"
#include "model.h"
#include "net.h"


/* SET UP REQUIRED RECORD SIZES.  Doesn't set the sizes of the variable-sized
   'S' and 'W' records. */

void net_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{
  logg->req_size['A'] = sizeof (net_arch);
  logg->req_size['F'] = sizeof (net_flags);
  logg->req_size['M'] = sizeof (model_specification);
  logg->req_size['V'] = sizeof (model_survival);
  logg->req_size['P'] = sizeof (net_priors);
  logg->req_size['D'] = sizeof (data_specifications);
}


/* REPORT ERROR IF SPECS FOR NET ARCHITECTURE, DATA MODEL, ETC. ARE MISSING. */

void net_check_specs_present
( net_arch *a,			/* Architecture, or null */
  net_priors *p,		/* Priors, or null */
  int need_model,		/* Must model be present? */
  model_specification *m,	/* Model, or null */
  model_survival *v		/* Survival model parameters, or null */
)
{
  if (a==0)
  { fprintf(stderr,"No architecture specification in log file\n");
    exit(1);
  }

  if (p==0)
  { fprintf(stderr,"No prior specification in log file\n");
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
