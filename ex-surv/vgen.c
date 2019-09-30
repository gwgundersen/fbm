/* VGEN.C - Generate survival data with hazard that varies in time. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 1000
#define N_censored 0
#define Censor_point 0.5

main()
{
  double x, h, u, t;
  int i;

  for (i = 0; i<N_cases; i++)
  {
    x = rand_gaussian();
    h = 1+(x-1)*(x-1);
    u = rand_uniopen();
    t = sqrt(-log(u)*2/h);

    if (i<N_censored && t>Censor_point) t = -Censor_point;

    printf(" %+8.5f %+10.5f\n",x,t);
  }
}
