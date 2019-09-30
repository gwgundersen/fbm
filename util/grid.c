/* GRID.C - Output grid of data points. */

/* Copyright (c) 1995-2004 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided 
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from Radford Neal's
 * web page, and note is made of any changes made to the programs.  The
 * programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
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

  exit(0);
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
    { sprintf(s,"%s %+.6e",prefix,val);
      grid(r+1,n-1,s);
    }
    else
    { printf("%s %+.6e\n",prefix,val);
    }
  }
}


/* PRINT USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf(stderr,"Usage: grid { low:high%%mod }\n");
  exit(1);
}
