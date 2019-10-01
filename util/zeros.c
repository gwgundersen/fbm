/* ZEROS.C - Create file of zeros. */

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

static void usage (void);

/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{ 
  int nlines, ncolumns;
  char junk;
  int i, j;

  if (argc!=3 || sscanf(argv[1],"%d%c",&nlines,&junk)!=1 || nlines<0
              || sscanf(argv[2],"%d%c",&ncolumns,&junk)!=1 || ncolumns<0) 
  { usage();
  }

  for (i = 0; i<nlines; i++)
  { for (j = 0; j<ncolumns; j++)
    { printf(" 0");
    }
    printf("\n");
  }

  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf(stderr,"Usage: zeros lines columns\n");
  exit(1);
}
