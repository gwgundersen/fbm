/* Generate data for the linear regression example. */

#include <stdio.h>
#include <math.h>
#include "rand.h"

#define N 100

main()
{
  double i0, i1, t;
  int i;

  for (i = 0; i<N; i++)
  { 
    i0 = rand_gaussian();
    i1 = i0 + 0.1*rand_gaussian();

    t = 0.5 + 2.5*i0 - 0.5*i1 + 0.1*rand_gaussian();

    printf("%10.5f %10.5f %10.5f\n",i0,i1,t);
  }

  exit(0);
}
