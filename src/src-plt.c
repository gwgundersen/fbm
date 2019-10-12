/* SRC-PLT.C - Procedures to plot data from source location model. */

/* Copyright (c) 2007 by Radford M. Neal 
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
#include "quantities.h"
#include "mc.h"
#include "src.h"


/* WHAT PROGRAM WE'RE PART OF. */

extern enum { PLT, TBL, HIST } program_type;


/* LOCAL PROCEDURES. */

extern void src_initialize (log_gobbled *);
extern void src_available  (quantities_described, log_gobbled *);
extern void src_evaluate   (quantities_described, quantities_held *,
                            log_gobbled *);
extern void src_mc_cleanup (void);

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
{ src_initialize,
  mc_initialize,
  0 
};

void (*quant_app_available[]) (quantities_described, log_gobbled *) = 
{ src_available,
  mc_available, 
  0 
};

void (*quant_app_evaluate[]) (quantities_described, quantities_held *,
                              log_gobbled *) =
{ src_evaluate,
  mc_evaluate,
  0
};

void (*quant_app_cleanup[]) (void) = 
{ src_mc_cleanup,
  0
};


/* DISPLAY USAGE MESSAGE AND EXIT. */

void plt_usage
( char *str
)
{ 
  fprintf(stderr,"Usage: src-%s\n",str);
}
