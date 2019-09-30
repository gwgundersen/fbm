/* DFT-PRINT.C - Procedures for printing parameters, etc. of diffusion tree. */

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
#include "dft.h"


/* PRINT HYPERPARAMETER VALUES. */

void dft_print_hypers
( dft_spec *dft,		/* Diffusion tree model specification */
  model_specification *m,	/* Specification for data model */
  dft_hypers *h			/* Structure containing hyperparameter values */
)
{
  int i, j;

  if (m && m->type=='R')
  { 
    printf("\nNOISE HYPERPARAMETERS\n\n");

    printf("%9.3f: ",h->noise_cm);

    j = dft->N_trees * dft->N_targets;

    for (i = 0; i<dft->N_targets; i++)
    { if (i>0 && i%5==0)
      { printf("\n           ");
      }
      printf(" %8.3f",h->SD[j++]);
    }

    printf("\n\n");
  }
}


/* PRINT PARAMETERS FOR TREES. */

void dft_print_params
( dft_spec *dft,		/* Diffusion tree model specification */
  dft_hypers *h			/* Structure containing hyperparameter values */
)
{
  int dt, i, j;

  printf("\n");

  j = 0;

  for (dt = 0; dt<dft->N_trees; dt++)
  { 
    printf("PARAMETERS OF TREE %d\n\n",dt+1);

    printf("Standard deviation parameters for diffusion process\n\n");

    printf("%9.3f: ",h->d_SD_cm[dt]);

    for (i = 0; i<dft->N_targets; i++)
    { if (i>0 && i%5==0)
      { printf("\n           ");
      }
      printf(" %8.3f",h->SD[j++]);
    }

    printf("\n\n");

    printf("Divergence function parameters:");
    if (h->c0[dt]==0) printf(" -"); else printf(" %.4f",h->c0[dt]);
    if (h->c1[dt]==0) printf(" -"); else printf(" %.4f",h->c1[dt]);
    if (h->c2[dt]==0) printf(" -"); else printf(" %.4f",h->c2[dt]);
    printf("\n\n");
  }
}


/* PRINT LATENT VECTORS FOR A SET OF CASES. */

void dft_print_latent
( dft_spec *dft,		/* Diffusion tree model specification */
  double *latent,		/* Latent vectors */
  int N_cases,			/* Number of cases */
  int bare			/* Print with no headings, more precision? */
)
{ 
  int i, j, k;

  if (!bare) printf("\nLATENT VECTORS\n\n");

  j = 0;

  for (i = 0; i<N_cases; i++)
  { printf ("%5d%c", i+1, bare ? ' ' : ':');
    for (k = 0; k<dft->N_targets; k++)
    { if (!bare && k!=0 && k%10==0) printf("\n      ");
      printf (bare ? " %+.8e" : " %+6.2f", latent[j++]);
    }
    printf("\n");
  }

  if (!bare) printf("\n");
}


/* PRINT TREES.  Show parents and divergence times. */

void dft_print_trees
( dft_spec *dft,		/* Diffusion tree model specification */
  dft_state st,			/* State variables */
  int N_cases,			/* Number of cases */
  int bare			/* Print without headings, and with tree nos? */
)
{
  int dt, i, j;

  j = 0;

  for (dt = 0; dt<dft->N_trees; dt++)
  {
    if (!bare)
    { printf("\nSTRUCTURE OF TREE %d\n",dt+1);
      printf("\n  Node  Parent  Div.Time\n\n");
    }
    
    for (i = -(N_cases-1); i<=N_cases; i++)
    { if (i==0) continue;
      if (bare) printf("%2d ",dt+1);
      printf ("%6d", i);
      printf ("  %6d", st[dt].parents[i]);
      printf ("  %8.6f\n", i<0 ? st[dt].divt[-i] : 1.0);
    }

  }

  if (!bare) printf("\n");
}


/* DISPLAY TREES IN GRAPHICAL FORM. */

static void graph
( dft_spec *dft,
  dft_state st,
  double time,
  int dt,
  int x,
  int s,
  int w
)
{ 
  int i, n0, n1, ca, cb;

  if (x>0)
  { printf(" %d",x);
  }
  else
  { 
    n0 = npts (st[dt].nodes[-x], 0);
    n1 = npts (st[dt].nodes[-x], 1);
    ca = chld (st[dt].nodes[-x], 0);
    cb = chld (st[dt].nodes[-x], 1);
    if (n0<n1 || n0==n1 && ca>cb)
    { ca = chld (st[dt].nodes[-x], 1);
      cb = chld (st[dt].nodes[-x], 0);
    }

    if (st[dt].divt[-x]<=time)
    {
      printf("%*d|",w,totpts(st[dt].nodes[-x]));

      graph (dft, st, time, dt, ca, s+1, w);
      printf("\n");
      for (i = 0; i<s+1; i++) printf("%*s|",w,"");
      if (cb<0)
      { printf("\n");
        for (i = 0; i<s+1; i++) printf("%*s|",w,"");
      }
      graph (dft, st, time, dt, cb, s+1, w);
    }
    else
    { graph (dft, st, time, dt, ca, s+1, w);
      graph (dft, st, time, dt, cb, s+1, w);
    }
  }
}

void dft_graph_trees
( dft_spec *dft,		/* Diffusion tree model specification */
  dft_state st,			/* State variables */
  int N_cases,			/* Number of cases */
  double time,			/* Time to go to */
  int bare			/* Print without headings? */
)
{
  int dt, root, w;

  w = 2 + (int) (log(N_cases+0.1) / log(10.0));

  for (dt = 0; dt<dft->N_trees; dt++)
  {
    if (!bare)
    { if (time<1)
      { printf(
         "\nGRAPH OF STRUCTURE OF TREE %d, FOR DIVERGENCE UP TO TIME %.5f\n\n",
         dt+1, time);
      }
      else
      { printf("\nGRAPH OF STRUCTURE OF TREE %d\n\n",dt+1);
      }
    }

    for (root = 1; st[dt].parents[root]!=0; root = st[dt].parents[root]) ;

    graph(dft,st,time,dt,root,0,w);
    printf("\n");
  }

  if (!bare) printf("\n");
}


/* PRINT NON-TERMINAL NODES.  Show children, number of terminal descendents, 
   divergence time, maybe location. */

void dft_print_nodes
( dft_spec *dft,		/* Diffusion tree model specification */
  dft_state st,			/* State variables */
  int N_cases,			/* Number of cases */
  int include_locations,	/* Print locations (if present)? */
  int bare			/* Print without headings, with tree nos, and
				   higher precision? */
)
{ 
  int dt, i, j;

  for (dt = 0; dt<dft->N_trees; dt++)
  {
    if (!bare)
    { printf("\nNON-TERMINAL NODES IN TREE %d\n",dt+1);
      printf("\n  Node  Children  Points  Div.Time");
      if (include_locations && st[dt].locations) printf("  Location");
      printf("\n\n");
    }

    for (i = 1; i<=N_cases-1; i++)
    { if (bare) printf("%2d ",dt+1);
      printf ("%6d %4d %4d  %6d", -i, 
              chld(st[dt].nodes[i],0), chld(st[dt].nodes[i],1),
              totpts(st[dt].nodes[i]));
      printf ("  %8.6f", st[dt].divt[i]);
      if (include_locations && st[dt].locations)
      { for (j = 0; j<dft->N_targets; j++)
        { printf (bare ? " %+.8e" : " %+6.2f", 
                  st[dt].locations[(i-1)*dft->N_targets+j]);
        }
      }
      printf("\n");
    }

  }

  if (!bare) printf("\n");
}
