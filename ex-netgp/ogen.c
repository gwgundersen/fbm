/* OGEN.C - Generate cases for a simple regression problem, with outliers. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 200

main()
{
  double x, y;
  int i;

  for (i = 0; i<N_cases; i++)
  {
    x = rand_gaussian();
    y = 0.3 + 0.4*x + 0.5*sin(2.7*x) + 1.1/(1+x*x);
    if (rand_uniform()<0.05)
    { y += rand_gaussian();
    }
    else
    { y += 0.1*rand_gaussian();
    }

    printf(" %+8.5f %+10.5f\n",x,y);
  }
}
