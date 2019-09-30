/* NUMIN-TEST.C - Program to test numeric input module. */

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

#include <stdio.h>

#include "numin.h"

main
( int argc,
  char **argv
)
{ 
  numin_source ns;
  double r[Max_items];
  int n, c, i, j;

  numin_spec(&ns,"data@1,0",1);

  do
  {
    if (argc<3 || (n = atoi(argv[2]))<=0) 
    { fprintf(stderr,"Usage: numin-test { specification n-items }\n");
      exit(1);
    }

    printf("\nSpecification: %s  Number of items: %d\n\n",argv[1],n);

    numin_spec(&ns,argv[1],n);
    c = numin_start(&ns);
    printf("  Number of lines: %d\n\n",c);
    for (i = 0; i<c; i++)
    { numin_read(&ns,r);
      for (j = 0; j<n; j++) printf("  %10.3le",r[j]);
      printf("\n");
    }
  
    numin_close(&ns);

    argv += 2;
    argc -= 2;

  } while (argc>1);

  printf("\n");

  exit(0);
}
