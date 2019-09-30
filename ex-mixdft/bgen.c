/* BGEN.C - Generate cases for binary mixture problem. */

#include <stdio.h>
#include <math.h>

#include "rand.h"

#define N_cases 1000

#define N_targets 10
#define N_comp 4

double M[N_comp][N_targets] = 
{ 0.1, 0.2, 0.2, 0.2, 0.2, 0.8, 0.8, 0.8, 0.8, 0.7,
  0.1, 0.8, 0.8, 0.8, 0.8, 0.2, 0.2, 0.2, 0.2, 0.7,
  0.9, 0.2, 0.2, 0.8, 0.8, 0.2, 0.2, 0.8, 0.8, 0.7,
  0.9, 0.8, 0.8, 0.2, 0.2, 0.8, 0.8, 0.2, 0.2, 0.7
};

main()
{
  int i, j;
  double *m;

  for (i = 0; i<N_cases; i++)
  {
    m = M[rand_int(N_comp)];

    for (j = 0; j<N_targets; j++)
    { 
      printf(" %d", rand_uniform()<m[j]);
    }

    printf("\n");
  }
}
