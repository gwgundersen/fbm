/* PHI.C - Gaussian density, CDF, and inverse CDF functions. */

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

#include "phi.h"


#define sqrt2   1.4142135623730950488
#define pi2     6.2831853071795864770 
#define sqrt2pi 2.5066282746310005024

double erfc (double);


/* GAUSSIAN DENSITY FUNCTION. */

double phi (double x)
{
  return exp(-x*x/2)/sqrt2pi;
}


/* GAUSSIAN CUMULATIVE DISTRIBUTION FUNCTION. */

double Phi (double x)
{
  return x<0 ? 0.5*erfc(-x/sqrt2) : 1 - 0.5*erfc(x/sqrt2);
}


/* INVERSE OF GAUSSIAN CUMULATIVE DISTRIBUTION FUNCTION.  See R. Thisted,
   Elements of Statistical Computing, page 332. */

double Phi_inverse (double p)
{
  double ix, x, c, d, v, t;
  int i;

  if (p<0.5) 
  { return -Phi_inverse(1-p);
  }

  if (p>0.99999999)
  { return 5.612002;
  }

  if (p>0.9)
  { v = -2*log(1-p);
    x = log(pi2*v);
    t = x/v + (2-x)/(v*v) + (-14+6*x-x*x)/(2*v*v*v);
    x = sqrt(v*(1-t)) * 1.01;
  }
  else
  { x = -log(1-p) - 0.6;
  }

  for (i = 0; i<20; i++)
  {
    c = Phi(x);
    d = phi(x);

    ix = (p-c) / (d - 0.5*(p-c)*x);
    x = x + ix;

    if (ix<1e-14 && ix>-1e-14) break;
  }

  if (x<-5.612002 || x>5.612002) 
  { abort();
  }

  return x;
}


/* TEST PROGRAM. */

#if 0

int main
( int argc,
  char **argv
)
{  
  double atof();
  double x;

  if (argc!=3) abort();

  x = atof(argv[1]);

  switch (*argv[2])
  {
    case 'd': printf("%.20f\n",phi(x)); break;
    case 'c': printf("%.20f\n",Phi(x)); break;
    case 'i': printf("%.20f\n",Phi_inverse(x)); break;
  }

  return 0;
}

#endif
