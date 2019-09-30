/* PHI.C - Gaussian density, CDF, and inverse CDF functions. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


#define sqrt2   1.4142135623730950488
#define pi2     6.2831853071795864770 
#define sqrt2pi 2.5066282746310005024

double erfc (double);


/* GAUSIAN DENSITY FUNCTION. */

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
  double x, c, d;
  int i;

  if (p<0.5) 
  { return -Phi_inverse(1-p);
  }

  if (p>0.9)
  { double v, t, x;
    v = -2*log(1-p);
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

    x = x + (p-c) / (d - 0.5*(p-c)*x);
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
