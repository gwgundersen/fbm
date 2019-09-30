/* DFT-DIV.C - Divergence functions and procedures. */

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
#include "rand.h"
#include "dft.h"


/* COMPUTE DIVERGENCE FUNCTION.  The value for a time of exactly t=1 is zero
   if divergence before then should have occurred, and one if the cumulative
   divergence to this point is finite. */

double dft_div 
( dft_hypers *h,		/* Hyperparameters, including divergence func.*/
  int dt,			/* Index of tree (from 0) */
  double t			/* Time argument, non-negative, <= 1 */
)
{ if (t>1) abort();
  if (t==1) 
  { return h->c1[dt]==0 && h->c2[dt]==0 ? 1 : 0;
  }
  return h->c0[dt] + h->c1[dt] / (1-t) + h->c2[dt] / ((1-t)*(1-t));
}


/* COMPUTE CUMULATIVE DIVERGENCE FUNCTION. */

double dft_cdiv 
( dft_hypers *h,		/* Hyperparameters, including divergence func.*/
  int dt,			/* Index of tree (from 0) */
  double t			/* Time argument, non-negative, <= 1 */
)
{ if (t>1) abort();
  return h->c0[dt] * t - (h->c1[dt]>0 ? h->c1[dt] * log(1-t) : 0)
                       + (h->c2[dt]>0 ? h->c2[dt] * (1/(1-t) - 1) : 0);
}


/* COMPUTE INVERSE CUMULATIVE DIVERGENCE FUNCTION.  Some special cases are
   done quickly, with the general case done using bisection.  An argument
   of INFINITY will result in a value of exactly 1. */

#define dft_icdiv_precision 40	/* Bits of precision desired */

double dft_icdiv 
( dft_hypers *h,		/* Hyperparameters, including divergence func.*/
  int dt,			/* Index of tree (from 0) */
  double e			/* Argument, non-negative */
)
{ 
  double low, high, mid;
  double c0, c1, c2, q;
  int i;

  if (e==INFINITY) return 1;

  c0 = h->c0[dt];
  c1 = h->c1[dt];
  c2 = h->c2[dt];

  if (c0==0 && c2==0)
  { 
    return 1 - exp (-e/c1);
  }

  else if (c0==0 && c1==0)
  { 
    return e / (e+c2);
  }

  else if (c1==0 && c2==0)
  { 
    return e>c0 ? 1 : e/c0;
  }
 
  else if (c1==0)
  { 
    q = c0 + c2 + e;
    return (q - sqrt(q*q-4*c0*e)) / (2*c0);
  }
 
  else /* general case */
  {
    low = 0; high = 1;
  
    for (i = 0; i<dft_icdiv_precision; i++)
    { mid = (low+high)/2;
      if (dft_cdiv(h,dt,mid)<e) 
      { low = mid;
      }
      else
      { high = mid;
      }
    }  
  
    return (low+high)/2;
  }
}


/* GENERATE PATH IN TREE TO NEW POINT.  Simulates generation of a new point
   to add to a tree.  The index of the node at the end of the segment from
   which divergence occurs is stored in "end"; the time of divergence is stored
   in "divt".  A divergence time of exactly 1 is avoided if the divergence
   function should diverge before then.  If the final argument is non-null,
   it gives the way to go at each branch point; if it is null, the way to
   go is chosen randomly, with the right probabilities. */

void dft_gen_path
( dft_hypers *hyp,	/* Hyperparameters for trees */
  dft_state st,		/* Pointers into state */
  int dt,		/* Number of tree (from 0) */
  int first,		/* First node in tree */
  int *end,		/* Place to store index of node at end of segment */
  double *divt,		/* Place to store divergence time */
  int *way		/* Way to go (0/1) at each branch, or null */
)
{
  double e, f;
  int x, c;

  do /* loop to avoid returning divergence time of 1 if not allowed */
  {
    x = first;
    e = 0;
   
    for (;;)
    { 
      if (x>0)
      { e += rand_exp();
        break;
      }
  
      c = totpts(st[dt].nodes[-x]);
  
      e += c*rand_exp();
  
      f = dft_cdiv(hyp,dt,st[dt].divt[-x]);
  
      if (e<=f) break;
     
      e = f;

      if (way!=0)
      { x = chld (st[dt].nodes[-x], way[-x]);
      }
      else
      { x = c*rand_uniform() < npts(st[dt].nodes[-x],0)
             ? chld(st[dt].nodes[-x],0) : chld(st[dt].nodes[-x],1);
      }
    }
  
    *end = x;

    *divt = dft_icdiv (hyp, dt, e);
    if (x<0 && *divt>st[dt].divt[-x]) /* in case roundoff has lead to this */
    { *divt = st[dt].divt[-x];
    }

  } while (*divt>=1 && (hyp->c1[dt]>0 || hyp->c2[dt]>0));
}
