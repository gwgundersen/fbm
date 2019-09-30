/* RGEN.C - Generate cases for a problem with a binary response. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 1000

main()
{
  double x1, x2, x3, x4, d;
  int i, c;

  for (i = 0; i<N_cases; i++)
  {
    x1 = rand_uniopen();
    x2 = rand_uniopen();
    x3 = rand_uniopen();
    x4 = rand_uniopen();

    d = sqrt ((x1-0.4)*(x1-0.4) + (x2-0.5)*(x2-0.5));

    c = d<0.35 ? 0 : 0.8*x1+1.8*x2<0.6 ? 1 : 2;

    printf(" %8.5f %8.5f %8.5f %8.5f %d\n",
      x1+0.1*rand_gaussian(), x2+0.1*rand_gaussian(), 
      x3+0.1*rand_gaussian(), x4+0.1*rand_gaussian(), 
      c);
  }
}
