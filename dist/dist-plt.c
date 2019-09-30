/* DIST-PLT.C - Procedures used to plot data for specified distribution. */

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
#include "quantities.h"
#include "mc.h"
#include "dist.h"


extern void dist_initialize (log_gobbled *);
extern void dist_available  (quantities_described, log_gobbled *);
extern void dist_evaluate   (quantities_described, quantities_held *,
                            log_gobbled *);

extern void mc_initialize  (log_gobbled *);
extern void mc_available   (quantities_described, log_gobbled *);
extern void mc_evaluate    (quantities_described, quantities_held *,
                            log_gobbled *);

void (*quant_app_arguments[]) (char ***) =
{ 0 
};

void (*quant_app_record_sizes[]) (log_gobbled *) =
{ mc_record_sizes,
  0 
};

void (*quant_app_initialize[]) (log_gobbled *) =
{ dist_initialize,
  mc_initialize,
  0 
};

void (*quant_app_available[]) (quantities_described, log_gobbled *) = 
{ dist_available,
  mc_available, 
  0 
};

void (*quant_app_evaluate[]) (quantities_described, quantities_held *,
                              log_gobbled *) =
{ dist_evaluate,
  mc_evaluate,
  0
};

void (*quant_app_cleanup[]) (void) = 
{ 0
};


/* DISPLAY USAGE MESSAGE AND EXIT. */

void plt_usage
( char *str
)
{ 
  fprintf(stderr,"Usage: dist-%s\n",str);
}
