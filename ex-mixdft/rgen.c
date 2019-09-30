/* RGEN.C - Generate cases for a bivariate density estimation problem. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 1000

main()
{
  double s, x, y;
  int i, j;

  for (i = 0; i<N_cases; i++)
  {
    if (rand_uniform()<0.3)
    { s = rand_gaussian();
      x = s+1.7;
      y = 1.9*s;
      for (j = 0; j<20; j++) 
      { s = rand_gaussian();
        y += s*s;
      }
    }
    else
    { s = rand_gaussian();
      x = s-1.9;
      y = -1.7*s;
      for (j = 0; j<10; j++) 
      { s = 0.3 + rand_gaussian();
        y += s*s;
      }
    }

    printf(" %+8.5lf %+10.5lf\n",x,y);
  }
}
