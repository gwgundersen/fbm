/* GRID.C - Output grid of data points. */

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
#include <string.h>
#include <stdio.h>
#include <math.h>


#define Max_dim 10

typedef struct { double low, high, mod; } range;

static void grid (range *, int, char *);
static void usage (void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{ 
  range r[Max_dim];
  int i, n;
  char *p;

  if (argc<2) usage();

  n = argc-1;

  if (n>Max_dim)
  { fprintf(stderr,"Grid has too many dimensions (max %d)\n",Max_dim);
    exit(1);
  }

  for (i = 0; i<n; i++)
  { p = argv[i+1];
    r[i].low = atof(p);
    while (*p!=':' && *p!=0) p += 1;
    if (*p==0) usage();
    p += 1;
    r[i].high = atof(p);
    while (*p!='%' && *p!=0) p += 1;
    if (*p==0) usage();
    p += 1;
    r[i].mod = atof(p);
    if (r[i].high<=r[i].low || r[i].mod<=0) usage();
  }

  grid(r,n,"");
}


/* OUTPUT THE GRID. */

static void grid
( range *r,
  int n,
  char *prefix
)
{
  char s[1000];
  double val;
  int i;

  for (i = ceil((r->low-1e-6)/r->mod); i*r->mod<=r->high+1e-6; i++)
  {
    val = i*r->mod;

    if (n>1)
    { sprintf(s,"%s %+.6le",prefix,val);
      grid(r+1,n-1,s);
    }
    else
    { printf("%s %+.6le\n",prefix,val);
    }
  }
}


/* PRINT USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf(stderr,"Usage: grid { low:high%%mod }\n");
  exit(1);
}
