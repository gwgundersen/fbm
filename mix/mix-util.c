/* MIX-UTIL.C - Utility routines for mixture modeling programs. */

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
#include "data.h"
#include "prior.h"
#include "model.h"
#include "mix.h"


/* SET UP REQUIRED RECORD SIZES.  Doesn't set the sizes of the variable-sized
   records. */

void mix_record_sizes
( log_gobbled *logg     /* Structure to hold gobbled data */
)
{
  logg->req_size['P'] = sizeof (mix_spec);
  logg->req_size['S'] = sizeof (mix_hypers);
  logg->req_size['M'] = sizeof (model_specification);
  logg->req_size['V'] = sizeof (model_survival);
}


/* REPORT ERROR IF MIXTURE MODEL SPECS, DATA MODEL, ETC. ARE MISSING. */

void mix_check_specs_present
( mix_spec *mx,			/* Mixture model specifications, or null */
  int need_model,		/* Must model be present? */
  model_specification *m	/* Model, or null */
)
{
  if (mx==0)
  { fprintf(stderr,"No specification for mixture model in log file\n");
    exit(1);
  }

  if (m==0 && need_model)
  { fprintf(stderr,"No model specification in log file\n");
    exit(1);
  }
}


/* SET UP INITIAL VALUES FOR HYPERPARAMETERS.  Sets the hyperparameters
   to the location parameters from their priors. */

void mix_hyper_init
( mix_spec *mx,			/* Specification for mixture model */
  model_specification *m,	/* Specification for data model */
  mix_hypers *h			/* Hyperparameter structure to initialize */
)
{
  int t;

  h->con = mx->con_prior.width;

  h->SD_cm = prior_width_scaled (&mx->SD_prior, mx->N_targets);

  for (t = 0; t<mx->N_targets; t++) 
  { h->SD[t] = h->SD_cm;
  }

  for (t = 0; t<mx->N_targets; t++)
  { h->mean[t] = 0;
  }

  if (m!=0 && m->type=='R')
  {
    h->noise_cm = m->noise.width;

    for (t = 0; t<mx->N_targets; t++)
    { h->noise[t] = h->noise_cm;
    }
  }
}


/* FIND FREQUENCIES OF COMPONENTS IN TRAINING CASES.  Looks at the indicators
   for the N_train training cases in the 'ind' array (numbered starting at 0),
   and stores frequencies of various components in 'freq'.  The 'freq' array 
   must be at least N_active long, and there must be no indicators bigger than 
   N_active-1. */

void mix_freq
( short *ind,		/* Array of component indicators for training cases */
  int N_train,		/* Number of training cases */
  int *freq,		/* Array to store frequencies of components */
  int N_active		/* Number of active components */
)
{
  int i, x;

  for (x = 0; x<N_active; x++)
  { freq[x] = 0;
  }

  for (i = 0; i<N_train; i++)
  { 
    x = ind[i];

    if (x>=N_active)
    { fprintf (stderr, "mix_freq: Oversize component indicator (%d>=%d)!\n",
               x, N_active);
      exit(1);
    }

    freq[x] += 1;
  }

}


/* SORT COMPONENTS BY DECREASING FREQUENCY.  The array of indicators, the 
   array of frequencies, and the arrays of component offsets and noise standard
   deviations (if they exist) are all updated to reflect the new component 
   numbering. */

void mix_sort
( short *ind,		/* Array of component indicators for training cases */
  int N_train,		/* Number of training cases */
  int *freq,		/* Array of frequencies of components */
  int N_active,		/* Number of active components */
  double *offsets,	/* Array of offsets for components, or null */
  double *noise_SD,	/* Array of noise SD for components, or null */
  int N_targets 	/* Number of target variables */
)
{
  double of[Max_targets], ns[Max_targets];
  int x, y, i, f, t;

  for (x = 1; x<N_active; x++)
  {
    f = freq[x];

    if (f>freq[x-1])
    {
      if (offsets) 
      { for (t = 0; t<N_targets; t++)
        { of[t] = offsets[x*N_targets+t];
        }
      }

      if (noise_SD)    
      { for (t = 0; t<N_targets; t++)
        { ns[t] = noise_SD[x*N_targets+t];
        }
      }

      for (i = 0; i<N_train; i++)
      { if (ind[i]==x) ind[i] = -1;
      } 

      for (y = x; y>0 && f>freq[y-1]; y--)
      {
        for (i = 0; i<N_train; i++)
        { if (ind[i]==y-1) ind[i] = y;
        }

        freq[y] = freq[y-1];

        if (offsets) 
        { for (t = 0; t<N_targets; t++)
          { offsets[y*N_targets+t] = offsets[(y-1)*N_targets+t];
          }
        }

        if (noise_SD)
        { for (t = 0; t<N_targets; t++)
          { noise_SD[y*N_targets+t] = noise_SD[(y-1)*N_targets+t];
          }
        }
      }

      freq[y] = f;

      if (offsets) 
      { for (t = 0; t<N_targets; t++)
        { offsets[y*N_targets+t] = of[t];
        }
      }

      if (noise_SD)    
      { for (t = 0; t<N_targets; t++)
        { noise_SD[y*N_targets+t] = ns[t];
        }
      }

      for (i = 0; i<N_train; i++)
      { if (ind[i]==-1) ind[i] = y;
      } 
    }
  }
}
