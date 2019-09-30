/* CMED.C - Generate true median times of death, corresponding to cgen.c. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 500

main()
{
  double x, h, t;
  int i;

  for (i = 0; i<N_cases; i++)
  {
    x = rand_gaussian();
    h = 1+(x-1)*(x-1);
    t = rand_exp()/h; /* Done to keep random numbers in sync */

    t = -log(0.5)/h;

    printf(" %+8.5f %+10.5f\n",x,t);
  }
}
