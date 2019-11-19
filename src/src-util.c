/* SRC-UTIL.C - Utility routines for source location modeling programs. */

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

#include "log.h"
#include "rand.h"
#include "src.h"


/* FIND LENGTH OF SOURCE PARAMETERS RECORD. */

int src_params_length
( int highN		/* Maximum number of sources */
)
{ src_params *junk_params;
  return ((char*)(junk_params->src+(highN)) - (char*)junk_params);
}


/* SET UP REQUIRED RECORD SIZES. */

void src_record_sizes
( log_gobbled *logg     /* Structure to hold gobbled data */
)
{
  logg->req_size['S'] = sizeof (src_spec);
  logg->req_size['T'] = sizeof (det_spec);
  logg->req_size['F'] = sizeof (flow_spec);
  logg->req_size['r'] = sizeof (rand_state);
}


/* SET DEFAULT INITIAL VALUES FOR PARAMETERS. */

void src_default_parameters 
( src_spec *src,		/* Source specification */
  det_spec *det,		/* Detector specification */
  flow_spec *flow,		/* Flow specification */
  src_params *params		/* Place to store initial values */
)
{
  int i, j;

  params->N0 = src->lowN;
  params->U = flow->lowU;
  params->log_width = det->log_low_width;
  params->inv_df = det->inv_low_df;

  for (i = 0; i<src->highN; i++)
  { for (j = 0; j<3; j++)
    { params->src[i].coord[j] = (src->low[j]+src->high[j])/2;
    }
    params->src[i].Q = 0;
    params->src[i].start = 0;
    params->src[i].stop = src->max_stop>=1e30 && src->max_duration>=1e30
                            ? 1e30 : 0;
  }

}


/* REPORT ERROR IF ANY SPECIFICATIONS ARE MISSING. */

void src_check_specs_present
( src_spec *src,		/* Source specification, or 0 */
  det_spec *det,		/* Detector specification, or 0 */
  flow_spec *flow		/* Flow specification, or 0 */
)
{
  if (src==0)
  { fprintf(stderr,"No source specification in log file\n");
    exit(1);
  }

  if (det==0)
  { fprintf(stderr,"No detector specification in log file\n");
    exit(1);
  }

  if (flow==0)
  { fprintf(stderr,"No flow specification in log file\n");
    exit(1);
  }
}
