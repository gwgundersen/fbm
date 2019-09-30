/* RGEN.C - Generate cases for a simple regression problem. */

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
    y = 0.3 + 0.4*x + 0.5*sin(2.7*x) + 1.1/(1+x*x) + 0.1*rand_gaussian();

    printf(" %+8.5f %+10.5f\n",x,y);
  }
}
