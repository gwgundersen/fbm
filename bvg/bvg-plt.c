/* BVG-PLT.C - Procedures used to plot data on bivariate Gaussian simulation. */

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
#include "quantities.h"
#include "mc.h"
#include "bvg.h"


/* WHAT PROGRAM WE'RE PART OF. */

extern enum { PLT, TBL, HIST } program_type;


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* LOCAL DATA. */

#define N_segments 1000	/* Number of segments plotted for ellipse */

static int ellipse;	/* Should we plot an ellipse? */


/* LOCAL PROCEDURES. */

static void bvg_arguments  (char ***);
static void bvg_initialize (log_gobbled *);

extern void mc_initialize  (log_gobbled *);
extern void mc_available   (quantities_described, log_gobbled *);
extern void mc_evaluate    (quantities_described, quantities_held *,
                            log_gobbled *);

void (*quant_app_arguments[]) (char ***) =
{ bvg_arguments,
  0 
};

void (*quant_app_record_sizes[]) (log_gobbled *) =
{ mc_record_sizes,
  0 
};

void (*quant_app_initialize[]) (log_gobbled *) =
{ bvg_initialize,
  mc_initialize,
  0 
};

void (*quant_app_available[]) (quantities_described, log_gobbled *) = 
{ mc_available, 
  0 
};

void (*quant_app_evaluate[]) (quantities_described, quantities_held *,
                              log_gobbled *) =
{ mc_evaluate,
  0
};

void (*quant_app_cleanup[]) (void) = 
{ 0
};


/* LOOK AT APPLICATON-SPECIFIC ARGUMENTS.  Only valid argument is "ellipse",
   and it's valid only for the plot program. */

static void bvg_arguments
( char ***argvp
)
{
  if (program_type==PLT && **argvp!=0 && strcmp(**argvp,"ellipse")==0)
  { 
    ellipse = 1;
    *argvp += 1;
  }
}


/* INITIALIZE FOR LOG FILE.  Only thing it may do is plot the ellipse for
   the distribution specified by this log file. */

static void bvg_initialize
( log_gobbled *logg
)
{
  static int cont = 0;

  bvg_spec *bs;
  double a, b;
  double x, y;
  int i;

  if (ellipse)
  { 
    if (cont) printf("\n");

    bs = logg->data['B'];
 
    if (bs==0)
    { fprintf(stderr,"No specification for bivariate Gaussian in log file\n");
      exit(1);
    }

    a = sqrt(1+bs->corr);
    b = sqrt(1-bs->corr);

    for (i = 0; i<=N_segments; i++)
    {
      x = a * cos(2.0*M_PI*i/N_segments);
      y = b * sin(2.0*M_PI*i/N_segments);

      printf ("%.6e %.6e\n", bs->std1*(x-y)/sqrt(2.0), 
                             bs->std2*(x+y)/sqrt(2.0));
    }

    printf("\n");
    cont = 1;
  }
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

void plt_usage
( char *str
)
{ 
  if (program_type==PLT)
  { fprintf(stderr,"Usage: bvg-%s [ / \"ellipse\" ]\n",str);
  }
  else
  { fprintf(stderr,"Usage: bvg-%s\n",str);
  }
}
