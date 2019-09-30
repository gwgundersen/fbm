/* RGEN.C - Generate cases for a problem with a binary response. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 500

main()
{
  double x1, x2, d;
  int i, b;

  for (i = 0; i<N_cases; i++)
  {
    x1 = rand_gaussian();
    x2 = rand_gaussian();

    d = (x1-0.4)*(x1-0.4) + (x2+0.3)*(x2+0.3);

    b = rand_uniform()<exp(-d*d);

    printf(" %+8.5f %+8.5f %d\n",x1,x2,b);
  }
}
