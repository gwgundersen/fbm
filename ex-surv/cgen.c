/* CGEN.C - Generate survival data with hazard that is constant in time. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 700
#define N_censored 200
#define Censor_point 0.5

main()
{
  double x, h, t;
  int i;

  for (i = 0; i<N_cases; i++)
  {
    x = rand_gaussian();
    h = 1+(x-1)*(x-1);
    t = rand_exp()/h;

    if (i<N_censored && t>Censor_point) t = -Censor_point;

    printf(" %+8.5f %+10.5f\n",x,t);
  }
}
