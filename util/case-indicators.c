/* CASE-INDICATORS.C - Output indicator variables identifying cases. */

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


static void usage (void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{ 
  int ntrn, ntst;
  char junk;
  int i, j;

  if (argc<2 || sscanf(argv[1],"%d%c",&ntrn,&junk)!=1 || ntrn<=0) usage();

  ntst = 0;
  if (argc>2)
  { if (sscanf(argv[2],"%d%c",&ntst,&junk)!=1 || ntst<0) usage();
  }
  
  if (argc>3) usage();

  for (i = 0; i<ntrn; i++)
  { for (j = 0; j<ntrn+1; j++)
    { printf(" %d",i==j);
    }
    printf("\n");
  }

  for (i = 0; i<ntst; i++)
  { for (j = 0; j<ntrn+1; j++)
    { printf(" %d",j==ntrn);
    }
    printf("\n");
  }

  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf(stderr,"Usage: case-indicators n-train [ n-test ]\n");
  exit(1);
}
