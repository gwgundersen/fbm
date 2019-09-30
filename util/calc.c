/* CALC.C - Calculate the value of a mathematical formula. */

/* Copyright (c) 1998 by Radford M. Neal 
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
#include <stdlib.h>
#include <string.h>

#include "formula.h"


static void usage(void)
{ fprintf(stderr,"Usage: calc { var=formula } formula \n");
  exit(1);
}

main
( int argc,
  char **argv
)
{
  char res[100];
  int a, c, i;
  char *e;

  if (argc<2) usage();

  for (a = 1; a<argc-1; a++)
  { 
    if (!strchr(argv[a],'='))
    { usage();
    }

    e = formula_def(argv[a],&c,&i);
    
    formula_var[c][i] = formula(e,0,1,0);
    formula_var_exists[c][i] = 1;
  }

  if (strchr(argv[argc-1],'=')) usage();

  sprintf (res, "%lg", formula(argv[argc-1],0,1,0));
  if (strchr(res,'e'))
  { *strchr(res,'e') = 'E';
  }

  printf("%s\n",res);

  exit(0);
}
