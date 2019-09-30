/* COMBINE.C - Combine files into one, with lines in parallel. */

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

#define Maxfiles 20

main(argc,argv)
  int argc;
  char **argv;
{ 
  FILE *f[Maxfiles];
  int nfiles, at_eof;
  int i, c;

  if (argc<2)
  { fprintf(stderr,"Usage: combine file1 file2 ...\n");
    exit(1);
  }
 
  nfiles = argc-1;

  for (i = 0; i<nfiles; i++)
  { f[i] = fopen(argv[i+1],"r");
    if (f[i]==NULL)
    { fprintf(stderr,"Can't open file %s\n",argv[i+1]);
      exit(1);
    }
  }

  for (;;)
  { 
    at_eof = 0;

    for (i = 0; i<nfiles; i++)
    { 
      c = getc(f[i]);
      if (c==EOF && i>0 && !at_eof || at_eof && c!=EOF)
      { fprintf(stderr,"Files do not have same number of lines\n");
        exit(1);
      }
      if (c==EOF) 
      { at_eof = 1;
        continue;
      }

      if (i>0) printf(" ");

      while (c!=EOF && c!='\n')
      { printf("%c",c);
        c = getc(f[i]);
      }

    }

    if (at_eof)
    { exit(0);
    }
    else
    { printf("\n");
    }
  }
}
