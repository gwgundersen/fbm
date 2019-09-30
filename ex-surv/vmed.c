/* VMED.C - Generate true median times of death, corresponding to vgen.c. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 1000

main()
{
  double x, h, u, t;
  int i;

  for (i = 0; i<N_cases; i++)
  {
    x = rand_gaussian();
    h = 1+(x-1)*(x-1);
    u = rand_uniopen(); /* Done to keep random numbers in sync */

    t = sqrt(-log(0.5)*2/h);

    printf(" %+8.5f %+10.5f\n",x,t);
  }
}
