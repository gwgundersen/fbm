/* RAND-GEN.C - Program to genenerate random numbers. */

/* Copyright (c) 1995-2005 by Radford M. Neal 
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


/* Usage:

     rand-gen seed generator { parameters } sample-size

Using the seed given, generates a sample of the given size using the indicated
random generator from the 'rand' module. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "rand.h"

static void usage (void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  int seed, sample_size, np;
  char *generator;
  double p1, p2;

  char **ap;
  double x;
  int i, n;

  if (argc<4) usage();

  if ((seed = atoi(argv[1]))==0 && strcmp(argv[1],"0")!=0) usage();

  generator = argv[2];

  if      (strcmp(generator,"uniform")==0)  np = 0;
  else if (strcmp(generator,"uniopen")==0)  np = 0;
  else if (strcmp(generator,"int")==0)      np = 1;
  else if (strcmp(generator,"poisson")==0)  np = 1;
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

  if (*ap==0 || (sample_size = atoi(*ap++))<=0) usage();

  if (*ap!=0) usage();

  rand_seed(seed);

  for (n = 0; n<sample_size; n++)
  {
    if      (strcmp(generator,"uniform")==0)  x = rand_uniform();
    else if (strcmp(generator,"uniopen")==0)  x = rand_uniopen();
    else if (strcmp(generator,"int")==0)      x = rand_int((int)p1);
    else if (strcmp(generator,"poisson")==0)  x = rand_poisson(p1);
    else if (strcmp(generator,"gaussian")==0) x = rand_gaussian();
    else if (strcmp(generator,"exp")==0)      x = rand_exp();
    else if (strcmp(generator,"cauchy")==0)   x = rand_cauchy();
    else if (strcmp(generator,"gamma")==0)    x = rand_gamma(p1);
    else if (strcmp(generator,"beta")==0)     x = rand_beta(p1,p2);
    else abort();

    printf("%.10e\n",x);
  }

  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

static void usage (void)
{
  fprintf(stderr,
   "Usage: rand-gen seed generator { parameters } sample-size\n");

  exit(1);
}
