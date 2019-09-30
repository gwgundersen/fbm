/* CALC.C - Calculate the value of a mathematical formula. */

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
