/* MODEL.C - Procedures dealing with modeling target value. */

/* Copyright (c) 1995, 1996 by Radford M. Neal 
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
#include "rand.h"


/* RETURN NUMBER OF TARGET VALUES REQUIRED.  No model at all, indicated
   by a null pointer, is treated the same a the 'R' model. */

int model_targets
( model_specification *m,	/* Model specification */
  int N_outputs			/* Number of outputs of function used by model*/
)
{ 
  if (m==0 || m->type==0) return N_outputs;

  switch (m->type)
  { 
    case 'B': case 'R': return N_outputs;

    case 'C': case 'V': return 1;

    default:
    { fprintf(stderr,"Unknown data model: %c\n",m->type);
      exit(1);
    }
  }
}
