/* RAND-TEST.C - Program to test random number generators. */

/* Copyright (c) 1995 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "rand.h"


#define Max_bins 1000		/* Maximum number of histogram bins */

static void usage();


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  int seed, sample_size, bins, np;
  double low, high;
  char *generator;
  double p1, p2;

  double mean, variance;
  double tmean, tvariance;
  int undef_mean, undef_variance;

  int count[Max_bins];
  int under, over;

  char **ap;
  double x;
  int i, n;

  bins = 0;

  if (argc<5) usage();

  if ((seed = atoi(argv[1]))==0 && strcmp(argv[1],"0")!=0) usage();

  generator = argv[2];

  if      (strcmp(generator,"uniform")==0)  np = 0;
  else if (strcmp(generator,"uniopen")==0)  np = 0;
  else if (strcmp(generator,"int")==0)      np = 1;
  else if (strcmp(generator,"gaussian")==0) np = 0;
  else if (strcmp(generator,"exp")==0)      np = 0;
  else if (strcmp(generator,"cauchy")==0)   np = 0;
  else if (strcmp(generator,"gamma")==0)    np = 1;
  else if (strcmp(generator,"beta")==0)     np = 2;
  else
  { fprintf(stderr,"Unknown generator: %s\n",generator);
    exit(1);
  }

  ap = argv+3;

  if (np>0) 
  { if (*ap==0 || (p1 = atof(*ap++))<=0) usage();
  }
  if (np>1) 
  { if (*ap==0 || (p2 = atof(*ap++))<=0) usage();
  }

  if (*ap==0 || strcmp(*ap++,"/")!=0) usage();

  if (*ap==0 || (sample_size = atoi(*ap++))<=0) usage();

  if (*ap!=0)
  { low = atof(*ap++);
    if (*ap==0) usage();
    high = atof(*ap++);
    if (high<=low) usage();
    if (*ap==0 || (bins = atoi(*ap++))<=0) usage();
    if (bins>Max_bins) 
    { fprintf(stderr,"Too many histogram bins\n");
      exit(1);
    }
  }

  if (*ap!=0) usage();

  printf("\nTest of %s(",generator);
  if (np>0) printf("%.4lf",p1);
  if (np>1) printf(",%.4lf",p2);
  printf(") generator using sample of size %d with seed %d\n\n",
   sample_size, seed);

  undef_mean = undef_variance = 0;

  if (strcmp(generator,"uniform")==0)  
  { tmean = 0.5;
    tvariance = 1.0/12.0;
  }
  else if (strcmp(generator,"uniopen")==0)  
  { tmean = 0.5;
    tvariance = 1.0/12.0;
  }
  else if (strcmp(generator,"int")==0)      
  { tmean = (p1-1)/2;
    tvariance = p1*p1/3.0 - p1/2.0 + 1/6.0 - tmean*tmean;
  }
  else if (strcmp(generator,"gaussian")==0) 
  { tmean = 0;
    tvariance = 1;
  }
  else if (strcmp(generator,"exp")==0) 
  { tmean = 1;
    tvariance = 1;
  }
  else if (strcmp(generator,"cauchy")==0)
  { undef_mean = 1;
    undef_variance = 1;
  }
  else if (strcmp(generator,"gamma")==0)
  { tmean = p1;
    tvariance = p1;
  }
  else if (strcmp(generator,"beta")==0)     
  { tmean = p1 / (p1+p2);
    tvariance = (p1*p2) / ((p1+p2)*(p1+p2)*(p1+p2+1));
  }
  else 
  { abort();
  }

  mean = 0;
  variance = 0;

  if (bins>0)
  { for (i = 0; i<bins; i++) count[i] = 0;
    under = over = 0;
  }

  rand_seed(seed);

  for (n = 0; n<sample_size; n++)
  {
    if      (strcmp(generator,"uniform")==0)  x = rand_uniform();
    else if (strcmp(generator,"uniopen")==0)  x = rand_uniopen();
    else if (strcmp(generator,"int")==0)      x = rand_int((int)p1);
    else if (strcmp(generator,"gaussian")==0) x = rand_gaussian();
    else if (strcmp(generator,"exp")==0)      x = rand_exp();
    else if (strcmp(generator,"cauchy")==0)   x = rand_cauchy();
    else if (strcmp(generator,"gamma")==0)    x = rand_gamma(p1);
    else if (strcmp(generator,"beta")==0)     x = rand_beta(p1,p2);
    else abort();

    mean += x;
    variance += x*x;

    if (bins>0)
    { if (x<low) 
      { under += 1;
      }
      else
      { i = (int) ((x-low)/((high-low)/bins));
        if (i>=bins) 
        { over += 1;
        }
        else 
        { count[i] += 1;
        }
      }
    }
  }

  mean /= sample_size;
  variance /= sample_size;
  variance -= mean*mean;

  printf("Sample mean:     %.4lf",mean);
  if (undef_mean)
  { printf(" (true value: undefined)\n");
  }
  else
  { printf(" (true value: %.4lf)\n",tmean);
  }

  printf("Sample variance: %.4lf",variance);
  if (undef_variance)
  { printf(" (true value: undefined)\n");
  }
  else
  { printf(" (true value: %.4lf)\n",tvariance);
  }
  printf("\n");

  if (bins!=0)
  { printf("Histogram:\n");
    printf("                    under : %8d  %.5lf\n\n", 
      under, (double)under / sample_size);
    for (i = 0; i<bins; i++)
    { printf("  %10.4lf - %10.4lf : %8d  %.5lf\n", 
        i*(high-low)/bins + low, (i+1)*(high-low)/bins + low, 
        count[i], (double)count[i] / sample_size);
    }
    printf("\n                     over : %8d  %.5lf\n", 
      over, (double)over / sample_size);
    printf("\n");
  }  

  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

static void usage()
{
  fprintf(stderr,
   "Usage: rand-test seed generator { parameters } / sample-size [ low high bins ]\n");

  exit(1);
}
