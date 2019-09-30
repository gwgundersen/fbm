/* DIST-QUANTITIES.C - Module defining quantities for specified distributions.*/

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
#include "mc.h"
#include "quantities.h"
#include "formula.h"
#include "data.h"
#include "dist.h"
#include "dist-data.h"


/* GLOBAL VARIABLES */

static dist_spec *dst;


/* INITIALIZE AFTER FIRST RECORDS READ. */

void dist_initialize
( log_gobbled *logg
)
{ 
  /* Check that specification record is present. */

  dst = logg->data['d'];
  if (dst==0)
  { fprintf(stderr,"No distribution specification in log file\n");
    exit(1);
  }

  /* Look at formulas to create variables. */

  (void) formula (dst->energy, 1, 0, 0);
  if (dst->Bayesian)
  { (void) formula (dst->energy+strlen(dst->energy)+1, 1, 0, 0);
  }

  /* Read data, if this is a Bayesian posterior. */

  if (dst->Bayesian)
  {
    data_spec = logg->data['D'];

    if (data_spec!=0)
    { dist_data_free();   
      dist_data_read();
    }
  }
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE. */

void dist_available 
( quantities_described qd,
  log_gobbled *logg
)
{ 
  char letter;
  int low, high, mod;
  int v;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    low  = qd[v].low;
    high = qd[v].high;
    mod = qd[v].modifier;

    if (letter && qd[v].available==0 && (strchr("it",letter)==0 || mod!=-1))
    {
      if (letter=='i')
      { qd[v].available = data_spec && mod<data_spec->N_inputs 
                            && low!=-1 && low<N_train && high<N_train ? 1 : -1; 
        if (high==-1) qd[v].high = N_train-1;
      }
      
      else if (letter=='t')
      { qd[v].available = data_spec && mod<data_spec->N_targets 
                            && low!=-1 && low<N_train && high<N_train ? 1 : -1; 
        if (high==-1) qd[v].high = N_train-1;
      }
 
      else if (letter=='l')
      { qd[v].available = 
          data_spec && low<N_train && high<N_train && mod==-1 ? 1 : -1; 
        if (low!=-1 && high==-1) qd[v].high = N_train-1;
      }

      else if (letter=='P')
      { qd[v].available = low==-1 && mod==-1 ? 1 : -1;
      }

      else if (strchr(STATE_VARS,letter))
      { qd[v].available = low==-1 && mod<=9 
          && formula_var_exists[letter-'a'][mod<0?10:mod] ? 1 : -1;
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE. */

void dist_evaluate 
( quantities_described qd, 
  quantities_held *qh,
  log_gobbled *logg
)
{ 
  int mod, low, high;
  char letter;
  int v, i;

  if (logg->data['q']==0 || logg->index['q']!=logg->last_index)
  { fprintf(stderr,"  records missing\n"); return;
  }

  dist_unpack_vars(logg->data['q']);

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    low  = qd[v].low;
    high = qd[v].high;
    mod  = qd[v].modifier;

    if (letter && !qh->updated[v] && (strchr("it",letter)==0 || mod!=-1))
    {
      if (letter=='i')
      { for (i = low; i<=high; i++)
        { qh->value[v][i-low] = train_inputs[data_spec->N_inputs*i+mod];
        }
        qh->updated[v] = 1;
        break;
      }
      
      else if (letter=='t')
      { for (i = low; i<=high; i++)
        { qh->value[v][i-low] = train_targets[data_spec->N_targets*i+mod];
        }
        qh->updated[v] = 1;
        break;
      }

      else if (letter=='l')
      { if (low==-1)
        { *qh->value[v] = N_train==0 ? 0 
                                     : dist_total_likelihood(dst,1.0,0)/N_train;
        }
        else
        { for (i = low; i<=high; i++)
          { qh->value[v][i-low] = dist_likelihood(dst,i,1.0,0);
          }
        }
        qh->updated[v] = 1;
        break;
      }

      else if (letter=='P')
      { *qh->value[v] = - dist_prior(dst,0);
        qh->updated[v] = 1;
        break;
      }

      else if (strchr(STATE_VARS,letter))
      { *qh->value[v] = formula_var[letter-'a'][mod<0?10:mod];
        qh->updated[v] = 1;
      }
    }
  }
}
