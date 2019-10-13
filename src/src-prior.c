/* SRC-PRIOR.C - Routines relating to prior for source model. */

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


/* SAMPLE PARAMETERS FROM THE PRIOR. */

void src_prior_generate
( src_params *params,	/* Place to store parameters generated  */
  src_spec *src,	/* Specification of source model */
  det_spec *det,	/* Specification of detector noise model */
  flow_spec *flow	/* Specification of flow model */
)
{
  int i;

  if (strchr("tT",flow->type))
  { params->U = flow->lowU + (flow->highU-flow->lowU) * rand_uniopen();
  }
  else
  { abort();
  }

  params->log_width = det->log_low_width 
    + (det->log_high_width-det->log_low_width) * rand_uniopen();

  params->inv_df = det->inv_high_df 
    + (det->inv_low_df-det->inv_high_df) * rand_uniopen();

  params->N0 = src->lowN + (src->highN-src->lowN+1) * rand_uniopen();

  for (i = 0; i<src->highN; i++)
  { src_prior_one_src (params, src, det, flow, i);
  }
}


/* SAMPLE PARAMETERS FOR ONE SOURCE FROM THE PRIOR. */

void src_prior_one_src
( src_params *params,	/* Parameters to modify */
  src_spec *src,	/* Specification of source model */
  det_spec *det,	/* Specification of detector noise model */
  flow_spec *flow,	/* Specification of flow model */
  int i			/* Which source to sample for */
)
{
  double upper_time;
  int j;

  if (i>=src->highN) abort();

  for (j = 0; j<3; j++) 
  { params->src[i].coord[j] = 
      src->low[j] + (src->high[j]-src->low[j]) * rand_uniopen();
  }

  params->src[i].Q = pow(src->lowQ,src->powQ) + rand_uniopen() 
                      * (pow(src->highQ,src->powQ) - pow(src->lowQ,src->powQ));

  params->src[i].start = src->max_start * rand_uniopen();

  upper_time = params->src[i].start + src->max_duration;
  if (upper_time>src->max_stop)
  { upper_time = src->max_stop;
  }
  params->src[i].stop = upper_time>=1e30 ? 1e30 
    : params->src[i].start + (upper_time-params->src[i].start) * rand_uniopen();
}
