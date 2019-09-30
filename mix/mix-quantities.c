/* MIX-QUANTITIES.C - Module defining quantities for mixture models. */

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
#include "prior.h"
#include "model.h"
#include "data.h"
#include "mix.h"
#include "mix-data.h"


/* MIXTURE MODEL SPECIFICATIONS. */

static mix_spec *mx;
static model_specification *model;

static int Max_active;
static int *freq;

static int have_train_data;


/* INITIALIZE AFTER FIRST RECORDS READ. */

void mix_initialize
( log_gobbled *logg
)
{ 
  /* Check that required specification records are present. */

  mx     = logg->data['P'];
  model  = logg->data['M'];

  mix_check_specs_present(mx,0,model);

  /* Read training and test data, if present. */

  have_train_data = 0;
  data_spec = logg->data['D'];
  N_train = 0;

  if (data_spec!=0)
  {
    have_train_data = 1;

    mix_data_free ();   
    mix_data_read (1, 1, mx, model);
  }

  Max_active = N_train;
  if (mx->N_components!=0 && mx->N_components<Max_active)
  { Max_active = mx->N_components;
  }

  if (freq!=0) free(freq);

  freq = chk_alloc (N_train, sizeof *freq);
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE. */

void mix_available 
( quantities_described qd,
  log_gobbled *logg
)
{ 
  char model_type = model ? model->type : 0;

  char letter;
  int low, high, mod;
  int v;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    low  = qd[v].low;
    high = qd[v].high;
    mod = qd[v].modifier;

    if (letter && qd[v].available==0 && (letter!='t' || mod!=-1))
    {
      switch (letter)
      {
        case 't':
        { if (mod==-1) break;
          qd[v].available = -1;
          if (mod>=0 && mod<mx->N_targets 
           && low!=-1 && low<N_train && high<N_train)
          { if (high==-1) qd[v].high = N_train-1;
            qd[v].available = 1;
          }
          break;
        }

        case 'u':
        { qd[v].available = -1;
          if  (low!=-1 && low<mx->N_targets && high<mx->N_targets)
          { if (high==-1) qd[v].high = mx->N_targets-1; 
            qd[v].available = 1;
          }
          break;
        }

        case 'v': case 'V':
        { qd[v].available = -1;
          if  (mod==-1 && low<mx->N_targets && high<mx->N_targets)
          { if (low!=-1 && high==-1) qd[v].high = mx->N_targets-1; 
            qd[v].available = 1;
          }
          break;
        }

        case 'n': case 'N':
        { qd[v].available = -1;
          if (model && model->type=='R' && (mod==-1 && low==-1
               || low!=-1 && low<mx->N_targets && high<mx->N_targets))
          { if (low!=-1 && high==-1) qd[v].high = mx->N_targets-1; 
            qd[v].available = 1;
          }
          break;
        }

#       if 0
        case 'K':
        { qd[v].available = -1;
          if (mod==-1 && low==-1 && high==-1)
          { qd[v].available = 1;
          }
          break;
        }
#       endif

        case 'c': case 'C':
        { qd[v].available = -1;
          if (mod>0 && low==-1 && high==-1 && N_train>0)
          { qd[v].available = 1;
          }
          break;
        }

        case 'a':
        { qd[v].available = -1;
          if (mod<100 && low==-1 && high==-1 && N_train>0)
          { qd[v].available = 1;
          }
          break;
        }

        case 'o':
        { qd[v].available = -1;
          if (mod>=0 && mod<mx->N_targets 
           && low!=-1 && low<N_train && high<N_train)
          { if (high==-1) qd[v].high = N_train-1;
            qd[v].available = 1;
          }
          break;
        }

        case 'h':
        { qd[v].available = -1;
          if (mod==-1 && low!=-1 && low<N_train && high<N_train)
          { if (high==-1) qd[v].high = N_train-1;
            qd[v].available = 1;
          }
          break;
        }

        default:
        { break;
        }
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE. */

void mix_evaluate 
( quantities_described qd, 
  quantities_held *qh,
  log_gobbled *logg
)
{ 
  char model_type = model ? model->type : 0;

  int mod, low, high;
  double *inputs;
  double *targets;
  double *offsets;
  double *noise_SD;
  short *indicators;
  mix_hypers *hyp;
  char letter;
  int N_active;
  int v, i;

  if (logg->data['S']==0 || logg->index['S']!=logg->last_index)
  { fprintf(stderr,"  records missing\n"); return;
  }

  hyp = logg->data['S'];

  indicators = logg->data['I']!=0 && logg->index['I']==logg->last_index
                   ? (short*) logg->data['I'] : 0;

  offsets    = logg->data['O']!=0 && logg->index['O']==logg->last_index
                   ? (double*) logg->data['O'] : 0;

  noise_SD   = logg->data['N']!=0 && logg->index['N']==logg->last_index
                   ? (double*) logg->data['N'] : 0;

  N_active = 0;
  if (indicators)
  { for (i = 0; i<N_train; i++)
    { if (indicators[i]+1>N_active) N_active = indicators[i]+1;
    }
  }

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    low  = qd[v].low;
    high = qd[v].high;
    mod  = qd[v].modifier;

    if (letter && !qh->updated[v] && (letter!='t' || mod!=-1))
    {
      inputs = train_inputs;
      targets = train_targets;

      if (mod<0) mod = 0;

      switch (letter)
      { 
        case 't':
        { if (mod==-1) break;
          for (i = low; i<=high; i++)
          { qh->value[v][i-low] = train_targets[data_spec->N_targets*i+mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'u':
        { if (mod>0 && offsets==0) break;
          for (i = low; i<=high; i++)
          { qh->value[v][i-low] = mod<=0 || mod>N_active ? hyp->mean[i]
                                    : offsets[(mod-1)*mx->N_targets+i];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'v': case 'V':
        { if (low==-1)
          { *qh->value[v] = hyp->SD_cm;
            if (letter=='V')
            { *qh->value[v] *= *qh->value[v];
            }
          }
          else
          { for (i = low; i<=high; i++)
            { qh->value[v][i-low] = hyp->SD[i];
              if (letter=='V')
              { qh->value[v][i-low] *= qh->value[v][i-low];
              }
            } 
          }
          qh->updated[v] = 1;
          break;
        }

        case 'n': case 'N':
        { if (mod<=0)
          { if (low==-1)
            { *qh->value[v] = hyp->noise_cm;
              if (letter=='N')
              { *qh->value[v] *= *qh->value[v];
              }
            }
            else
            { for (i = low; i<=high; i++)
              { qh->value[v][i-low] = hyp->noise[i];
                if (letter=='N')
                { qh->value[v][i-low] *= qh->value[v][i-low];
                }
              }
            }
          }
          else
          { if (noise_SD==0) break;
            for (i = low; i<=high; i++)
            { qh->value[v][i-low] = mod>N_active ? hyp->noise[i]
                                       : noise_SD[(mod-1)*mx->N_targets+i];
              if (letter=='N')
              { qh->value[v][i-low] *= qh->value[v][i-low];
              }
            }
          }
          qh->updated[v] = 1;
          break;
        }

#       if 0
        case 'K':
        { *qh->value[v] = hyp->con;
          qh->updated[v] = 1;
          break;
        }
#       endif

        case 'c': case 'C':
        { if (indicators==0) break;
          mix_freq (indicators, N_train, freq, Max_active);
          *qh->value[v] = mod>Max_active ? 0 : freq[mod-1];
          if (letter=='C')
          { for (i = 1; i<mod && i<=Max_active; i++)
            { *qh->value[v] += freq[i-1];
            }
          }
          *qh->value[v] /= N_train;
          qh->updated[v] = 1;
          break;
        }

        case 'a':
        { if (indicators==0) break;
          if (mod<=0) 
          { *qh->value[v] = N_active;
          }
          else
          { int need;
            need = N_train - (int)(0.5+N_train*mod/100.0);
            mix_freq (indicators, N_train, freq, Max_active);
            for (i = 0; need>0; i++)
            { if (i>=Max_active) abort();
              need -= freq[i];
            }
            *qh->value[v] = i;
          }
          qh->updated[v] = 1;
          break;
        }

        case 'o':
        { if (indicators==0 || offsets==0) break;
          for (i = low; i<=high; i++)
          { qh->value[v][i-low] = offsets[indicators[i]*mx->N_targets+mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'h':
        { if (indicators==0) break;
          mix_freq(indicators,N_train,freq,N_active);
          for (i = low; i<=high; i++)
          { if (indicators[i]<0 || indicators[i]>=N_active) abort();
            qh->value[v][i-low] = freq[indicators[i]];
          }
          qh->updated[v] = 1;
          break;
        }
 
        default:
        { break;
        }
      }
    }
  }
}
