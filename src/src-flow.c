/* SRC-FLOW.C - Flow computations. */

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

#include "phi.h"
#include "log.h"
#include "src.h"


/* COMPUTE TOTAL FROM ALL SOURCES, FOR ANY FLOW TYPE KNOWN. */

double src_total
( src_spec *src,		/* Specification for source model */
  flow_spec *flow,		/* Specification of flow model */
  src_params *params,		/* Parameters of model */
  double rcoords[4]		/* Receptor coordinates + time */
)
{
  double v, c;
  int N, i;

  N = (int) params->N0;

  v = 0;

  for (i = 0; i<N; i++)
  { 
    if (flow->type=='t')
    { c = src_test_cstar (flow, params, params->src[i].coord, rcoords);
    }
    else if (flow->type=='T')
    { c = src_testss_cstar (flow, params, params->src[i].coord, 
                            params->src[i].start, params->src[i].stop, rcoords);
    }
    else
    { abort();
    }
    v += c * pow (params->src[i].Q, 1/src->powQ);
  }

  return v;
}


/* COMPUTE THE ADJUNCT FUNCTION FOR THE TEST FLOW MODEL. */

double src_test_cstar
( flow_spec *flow,		/* Specification of flow model */
  src_params *params,		/* Parameters of model */
  double scoords[3],		/* Source coordinates */
  double rcoords[3]		/* Receptor coordinates */
)
{
  double sigma_y, sigma_z;
  double xd, G, H;

  xd = scoords[0]-rcoords[0];
  if (xd<=0) 
  { return 0;
  }

  sigma_y = flow->ay * xd / sqrt (1+flow->by*xd);
  sigma_z = flow->az * xd / sqrt (1+flow->bz*xd);

  G = phi((scoords[1]-rcoords[1])/sigma_y) / sigma_y;
  H = (phi((scoords[2]-rcoords[2])/sigma_z) + 
       phi((scoords[2]+rcoords[2])/sigma_z)) / sigma_z;

  return (1/params->U) * G * H;
}


/* COMPUTE THE ADJUNCT FUNCTION FOR THE TEST-START-STOP FLOW MODEL. */

double src_testss_cstar
( flow_spec *flow,		/* Specification of flow model */
  src_params *params,		/* Parameters of model */
  double scoords[3],		/* Source coordinates */
  double start,			/* Source start time */
  double stop,			/* Source stop time */
  double rcoords[4]		/* Receptor coordinates + time */
)
{
  double sigma_y, sigma_z;
  double xd, t, G, H;

  xd = scoords[0]-rcoords[0];
  t = rcoords[3]-xd/params->U;

  if (xd<=0  || t<start || t>stop)
  { return 0;
  }

  sigma_y = flow->ay * xd / sqrt (1+flow->by*xd);
  sigma_z = flow->az * xd / sqrt (1+flow->bz*xd);

  G = phi((scoords[1]-rcoords[1])/sigma_y) / sigma_y;
  H = (phi((scoords[2]-rcoords[2])/sigma_z) + 
       phi((scoords[2]+rcoords[2])/sigma_z)) / sigma_z;

  return (1/params->U) * G * H;
}
