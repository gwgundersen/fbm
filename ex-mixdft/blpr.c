/* BLPR.C - Find log probabilities of cases w.r.t. the true mixture. */

#include <stdio.h>
#include <math.h>

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
  double pr[N_comp], prt;
  int j, k, b;

  for (;;)
  { 
    for (k = 0; k<N_comp; k++) 
    { pr[k] = 1;
    }

    for (j = 0; j<N_targets; j++)
    { if (scanf("%d",&b)!=1) 
      { if (j==0) 
        { exit(0);
        }
        else
        { abort();
        }
      }
      for (k = 0; k<N_comp; k++)
      { pr[k] *= b==1 ? M[k][j] : 1-M[k][j];
      }
    }

    prt = 0;
    for (k = 0; k<N_comp; k++)
    { prt += pr[k];
    }
    prt /= N_comp;

    printf("%.3f\n",log(prt));
  }
}
