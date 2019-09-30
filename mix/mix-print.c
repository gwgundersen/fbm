/* MIX-PRINT.C - Procedures for printing parameters of mixture model. */

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
#include "prior.h"
#include "model.h"
#include "data.h"
#include "mix.h"


/* PRINT HYPERPARAMETER VALUES. */

void mix_print_hypers
( mix_spec *mx,			/* Mixture model specification */
  model_specification *m,	/* Specification for data model */
  mix_hypers *h			/* Structure containing hyperparameter values */
)
{
  int i;

  printf("\nHYPERPARAMETERS\n");

# if 0

  printf("\nDirichlet concentration parameter: %.3f (unscaled)",h->con);

  if (mx->N_components>0 && mx->con_prior.scale) 
  { printf(", %.3f (scaled)",h->con/mx->N_components);
  }

  printf("\n");

# endif

  printf("\nStandard deviations for component offsets:\n");

  printf("\n%9.3f: ",h->SD_cm);

  for (i = 0; i<mx->N_targets; i++)
  { if (i>0 && i%5==0)
    { printf("\n           ");
    }
    printf(" %8.3f",h->SD[i]);
  }

  printf("\n");

  printf("\nMeans for component offsets:\n");

  for (i = 0; i<mx->N_targets; i++)
  { if (i%5==0)
    { printf("\n           ");
    }
    printf(" %+8.3f",h->mean[i]);
  }

  printf("\n");

  if (m && m->type=='R')
  { 
    printf("\nStandard deviations for Gaussian target distributions:\n");

    printf("\n%9.3f: ",h->noise_cm);

    for (i = 0; i<mx->N_targets; i++)
    { if (i>0 && i%5==0)
      { printf("\n           ");
      }
      printf(" %8.3f",h->noise[i]);
    }

    printf("\n");
  }
}


/* PRINT VALUES OF COMPONENT PARAMETERS (AND FREQUENCIES). */

void mix_print_params
( mix_spec *mx,			/* Mixture model specification */
  model_specification *m,	/* Specification of data model */
  int N_active,			/* Number of active components */
  int *freq,			/* Frequencies of components in dataset */
  double *offsets,		/* Offset parameters */
  double *noise_SD		/* Noise standard deviations, or null */
)
{
  int i, T, t;

  printf("\nPARAMETERS AND FREQUENCIES FOR COMPONENTS OF THE MIXTURE\n");

  T = 0;

  for (i = 0; i<N_active; i++)   
  { T += freq[i];
  }

  for (i = 0; i<N_active; i++)
  {
    printf ("\n%4d: %5.3f", i+1, (double)freq[i]/T);

    for (t = 0; t<mx->N_targets; t++)
    { if (t>0 && t%5==0) printf("\n           ");
      printf(" %+8.3f",offsets[i*mx->N_targets+t]);
    }

    printf ("\n");

    if (m->type=='R')
    { 
      for (t = 0; t<mx->N_targets; t++)
      { if (t%5==0) printf("\n           ");
        printf(" %8.3f",noise_SD[i*mx->N_targets+t]);
      }

      printf ("\n");
    }
  }
}


/* PRINT VALUES OF INDICATORS. */

void mix_print_indicators
( int N_cases,			/* Number of cases in indicator array */
  short *indicators		/* Array of indicators */
)
{
  int i;
 
  printf("\nINDICATORS OF WHICH COMPONENT IS ASSOCIATED WITH EACH CASE\n");

  for (i = 0; i<N_cases; i++)
  { 
    if (i%10==0) printf("\n%5d:",i);

    printf(" %5d",indicators[i]+1);
  }

  printf("\n");
}
