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
    case 'B': case 'N': case 'R': return N_outputs;

    case 'C': case 'V': return 1;

    default:
    { fprintf(stderr,"Unknown data model: %c\n",m->type);
      exit(1);
    }
  }
}
