/* SRC-QUANTITIES.C - Module defining quantities for source location model. */

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
#include "numin.h"
#include "data.h"
#include "src-data.h"


/* VARIABLES LOCAL TO THIS MODULE. */

static src_spec *src;
static det_spec *det;
static flow_spec *flow;


/* INITIALIZE AFTER FIRST RECORDS READ. */

void src_initialize
( log_gobbled *logg
)
{ 
  src = logg->data['S'];
  det = logg->data['T'];
  flow = logg->data['F'];

  src_check_specs_present (src, det, flow);
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE. */

void src_available 
( quantities_described qd,
  log_gobbled *logg
)
{ 
  char letter;
  int mod, low, high;
  int v;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    mod = qd[v].modifier;
    low = qd[v].low;
    high = qd[v].high;

    if (letter && qd[v].available==0)
    {
      if (letter=='N' && (mod==-1 || mod==0) && low==-1)
      { qd[v].available = 1;
      }
      else if (strchr("xyz",letter) && mod==-1 && low!=-1)
      { qd[v].available = high>=src->highN ? -1 : 1;
        if (qd[v].available && high==-1) qd[v].high = src->highN-1;
      }
      else if (letter=='Q' && (mod==-1 || mod==1) && low!=-1)
      { qd[v].available = high>=src->highN ? -1 : 1;
        if (qd[v].available && high==-1) qd[v].high = src->highN-1;
      }
      else if (letter=='R' && (mod==-1 || mod==1) && low!=-1)
      { qd[v].available = src->max_stop==1e30 && src->max_duration==1e30
                           || high>=src->highN ? -1 : 1;
        if (qd[v].available && high==-1) qd[v].high = src->highN-1;
      }
      else if (letter=='T' && (mod==0 || mod==1 || mod==2) && low!=-1)
      { qd[v].available = strchr(steady_state,flow->type) 
            || high>=src->highN 
            || mod>0 && src->max_stop==1e30 && src->max_duration==1e30
          ? -1 : 1;
        if (qd[v].available && high==-1) qd[v].high = src->highN-1;
      }
      else if (letter=='U'  && mod==-1 && low==-1 && flow->type=='t')
      { qd[v].available = 1;
      }
      else if (letter=='n' && mod==-1 && low==-1)
      { qd[v].available = 1;
      }
      else if (letter=='u' && mod==-1 && low==-1 && det->inv_low_df>0)
      { qd[v].available = 1;
      }
      else if (strchr("vo",letter) && mod==-1 && low!=-1)
      { qd[v].available = high>=N_train ? -1 : 1;
        if (qd[v].available && high==-1) qd[v].high = N_train;
      }
      else if (strchr("VO",letter) && mod==-1 && low!=-1)
      { qd[v].available = high>=N_test ? -1 : 1;
        if (qd[v].available && high==-1) qd[v].high = N_test;
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE. */

static int dbl_cmp (const void *a, const void *b) 
{ double d;
  d = *(double*)b - *(double*)a;
  return d==0 ? 0 : d<0 ? -1 : 1;
}

void src_evaluate 
( quantities_described qd, 
  quantities_held *qh,
  log_gobbled *logg
)
{ 
  src_params *params;
  int mod, low, high;
  char letter;
  int N, i, v;

  if (logg->data['q']==0 || logg->index['q']!=logg->last_index)
  { fprintf(stderr,"Parameter record missing\n"); 
    return;
  }

  params = logg->data['q'];

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    mod = qd[v].modifier;
    low = qd[v].low;
    high = qd[v].high;

    N = (int)params->N0;

    if (letter && !qh->updated[v])
    {
      if (letter=='N')
      { *qh->value[v] = params->N0;
        if (mod==-1) *qh->value[v] = (int)*qh->value[v];
        qh->updated[v] = 1;
      }

      else if (letter=='x' || letter=='y' || letter=='z')
      { 
        for (i = low; i<=high; i++)
        { qh->value[v][i-low] = params->src[i].coord[letter-'x'];
        }

        qh->updated[v] = 1;
      }

      else if (letter=='Q' || letter=='R')
      { 
        if (mod==1) 
        { double *Qvals;
          Qvals = N==0 ? 0 : chk_alloc (N, sizeof (double));
          for (i = 0; i<N; i++)
          { Qvals[i] = pow (params->src[i].Q, 1/src->powQ);
            if (letter=='R') 
            { Qvals[i] *= params->src[i].stop - params->src[i].start;
            }
          }
          qsort (Qvals, N, sizeof (double), dbl_cmp);
          for (i = low; i<=high; i++)
          { qh->value[v][i-low] = i<N ? Qvals[i] : 0;
          }
          free(Qvals);
        }
        else
        { for (i = low; i<=high; i++)
          { qh->value[v][i-low] = pow (params->src[i].Q, 1/src->powQ);
            if (letter=='R')
            { qh->value[v][i-low] *= params->src[i].stop - params->src[i].start;
            }
          }
        }

        qh->updated[v] = 1;
      }

      else if (letter=='T')
      { 
        for (i = low; i<=high; i++)
        { qh->value[v][i-low] = mod==0 ? params->src[i].start 
                              : mod==1 ? params->src[i].stop
                              : params->src[i].stop - params->src[i].start;
        }

        qh->updated[v] = 1;
      }

      else if (letter=='U')
      { 
        *qh->value[v] = params->U;
        qh->updated[v] = 1;
      }

      else if (letter=='n')
      { 
        *qh->value[v] = exp(params->log_width);
        qh->updated[v] = 1;
      }

      else if (letter=='u')
      { 
        *qh->value[v] = 1/params->inv_df;
        qh->updated[v] = 1;
      }

      else if (letter=='v' || letter=='V')
      {
        for (i = low; i<=high; i++)
        { qh->value[v][i-low] = letter=='m' ? train_targets[i]
                                            : test_targets[i];
        }

        qh->updated[v] = 1;
      }

      else if (letter=='o' || letter=='O')
      {
        for (i = low; i<=high; i++)
        { qh->value[v][i-low] = 
              letter=='o' ? src_total (src, flow, params, train_inputs+4*i)
                          : src_total (src, flow, params, test_inputs+4*i);
        }

        qh->updated[v] = 1;
      }
    }
  }
}
