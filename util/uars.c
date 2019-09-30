/* UARS.C - Adaptive rejection sampling for unimodal distributions on (0,1). */

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rand.h"
#include "uars.h"


#define Max_rectangles 100	/* Maximum no. of rectangles in approximation */


/* ADAPTIVE REJECTION SAMPLING FOR A UNIMODAL DISTRIBUTION ON (0,1).

   This procedure samples from an arbitrary unimodal distribution on the
   interval (0,1), using an adaptive rejection technique.

   The first argument, logp, is the function for evaluating the log density, 
   up to an arbitrary additive constant.  It takes two arguments: a real value
   in the interval (0,1), and a pointer to possible other parameters of the
   density.  The second argument is the location of the distribution's mode.
   The third argument is a pointer to the extra parameters to be passed to 
   logp each time it is called.

   The region under the curve of the unnormalized density is approximated
   by a union of rectangles, from which points are sampled.  Rejected points
   are used to improve the approximation, until Max_rectangles rectangles 
   have been created (at which point the approximation stays fixed).  A warning 
   is issued (one time only) when 10000 rejections are done for one sampling
   call.
*/

double uars
( double (*logp)(double,void*),	/* Log of unnormalized density function */
  double mode,			/* Location of mode */
  void *extra			/* Extra parameters passed to logp */
)
{ 
  double lowx[Max_rectangles+1], height[Max_rectangles], area[Max_rectangles];
  int nrect, nrej, i, j;
  double x, y, pr, lpm;

  static int warned = 0;

  if (mode>1 || mode<0) abort();

  /* Set up initial two rectangles (one on each side of the mode). */

  lpm = logp(mode,extra);

  lowx[0] = 0;
  height[0] = 1;
  area[0] = mode;

  lowx[1] = mode;
  height[1] = 1;
  area[1] = (1-mode);

  lowx[2] = 1;

  nrect = 2;

  /* Sample points from the rectangles until one is accepted. */

  nrej = 0;

  for (;;)
  {
    /* Sample x and y uniformly from the collection of rectangles. */

    i = rand_pickd(area,nrect);
    x = lowx[i] + (lowx[i+1] - lowx[i]) * rand_uniopen();
    y = height[i] * rand_uniform();

    /* See if x is acceptable, returning it if so. */

    pr = exp(logp(x,extra)-lpm);
    if (pr>height[i]*1.0001) abort();

    if (y<pr) 
    { return x;
    }
    else
    { nrej += 1;
    }

    /* Improve the approximation by rectangles using the rejected point,
       unless we've already reached the maximum number of rectangles. */

    if (nrect<Max_rectangles)
    { 
      for (j = nrect; j>i; j--)
      { lowx[j] = lowx[j-1];
        height[j] = height[j-1];
        area[j] = area[j-1];
      }

      nrect += 1;
      lowx[nrect] = 1;

      lowx[i+1] = x;
      if (x>mode)
      { height[i+1] = pr;
      }
      else
      { height[i] = pr;
      }
      area[i] = (lowx[i+1]-lowx[i]) * height[i];
      area[i+1] = (lowx[i+2]-lowx[i+1]) * height[i+1];
    }

    /* Warn (just once) if we've rejected 10000 times in one call. */

    if (nrej==10000 && !warned)
    { fprintf(stderr, "WARNING: More than 10000 rejections in uars\n");
      warned = 1;
    }
  }
}


/* TEST PROGRAM.  Samples from a Beta distribution.  */

#if 0

static double logp (double x, void *v) 
{ double *p = v;
  return (p[0]-1)*log(x) + (p[1]-1)*log(1-x); 
}

int main 
( int argc,
  char **argv
)
{
  double atof();
  double p[2];
  double mode;
  int n;
  
  if (argc!=4) abort();

  p[0] = atof(argv[1]);
  p[1] = atof(argv[2]);
  n = atoi(argv[3]);

  if (p[0]<1 || p[1]<1) abort();

  mode = p[0]==1 && p[1]==1 ? 0.5 : (p[0]-1)/(p[0]+p[1]-2);

  while (n>0)
  { printf ("%.6f\n", uars (logp, mode, p));
    n -= 1;
  }
}

#endif
