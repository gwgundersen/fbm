/* DEMO-PLT.C - Lists of procedures used in demo of plot program (and others).*/

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

#include "log.h"
#include "quantities.h"


/* LISTS OF PROCEDURES FOR HANDLING QUANTITIES. */

extern void demo_arguments    (char ***);
extern void demo_initialize   (log_gobbled *);
extern void demo_record_sizes (log_gobbled *);
extern void demo_available    (quantities_described, log_gobbled *);
extern void demo_evaluate     (quantities_described, quantities_held *,
                               log_gobbled *);

void (*quant_app_arguments[]) (char ***) =
{ demo_arguments,
  0 
};

void (*quant_app_record_sizes[]) (log_gobbled *) =
{ demo_record_sizes,
  0 
};

void (*quant_app_initialize[]) (log_gobbled *) =
{ demo_initialize,
  0 
};

void (*quant_app_available[]) (quantities_described, log_gobbled *) = 
{ demo_available, 
  0 
};

void (*quant_app_evaluate[]) (quantities_described, quantities_held *,
                              log_gobbled *) =
{ demo_evaluate,
  0
};

void (*quant_app_cleanup[]) (void) = 
{ 0
};


/* DISPLAY USAGE MESSAGE. */

void plt_usage
( char *str
)
{ 
  fprintf(stderr,"Usage: demo-%s [ / integer ]\n",str);
}
