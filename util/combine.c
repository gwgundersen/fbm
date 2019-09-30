/* COMBINE.C - Combine files into one, with lines in parallel. */

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

#define Maxfiles 100

main(argc,argv)
  int argc;
  char **argv;
{ 
  FILE *f[Maxfiles];
  double num[Maxfiles];
  int nfiles, at_eof;
  int any_files;
  char junk;
  int i, c;

  if (argc<2)
  { fprintf(stderr,"Usage: combine file1 file2|number2 ...\n");
    exit(1);
  }
 
  nfiles = argc-1;

  if (nfiles>Maxfiles)
  { fprintf(stderr,"Too many files/numbers\n");
    exit(1);
  }

  any_files = 0;

  for (i = 0; i<nfiles; i++)
  { if (i>0 && sscanf(argv[i+1],"%lf%c",&num[i],&junk)==1)
    { f[i] = NULL;
    }
    else
    { f[i] = fopen(argv[i+1],"r");
      if (f[i]==NULL)
      { fprintf(stderr,"Can't open file %s\n",argv[i+1]);
        exit(1);
      }
      any_files = 1;
    }
  }

  if (!any_files)
  { fprintf(stderr,"At least one argument must be a file\n");
    exit(1);
  }

  for (;;)
  { 
    at_eof = 0;

    for (i = 0; i<nfiles; i++)
    { 
      if (f[i]==NULL)
      { if (!at_eof) 
        { if (num[i]==(int)num[i])
          { printf(" %d",(int)num[i]);
          }
          else
          { printf(" %f",num[i]);
          }
        }
      }
      else
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
    }

    if (at_eof)
    { exit(0);
    }
    else
    { printf("\n");
    }
  }
}
