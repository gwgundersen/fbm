/* EXTRACT.C - Extract items at random from data file. */

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

#include "rand.h"


main(argc,argv)
  int argc;
  char **argv;
{ 
  char *full_file, *extract_file, *remainder_file;
  FILE *full, *extract, *remainder;
  double frac;
  int seed;
  int ch;

  int total, need, left, chosen;

  double atof();

  /* Look at arguments. */

  seed = 1;

  if (argc<5 || argc>6
   || argc>5 && (seed = atoi(argv[5]))<=0
   || (frac = atof(argv[1]))<=0 || frac>1.0)
  { fprintf(stderr,
      "Usage: extract frac full-file extract-file remainder-file [ seed ]\n");
    exit(-1);
  }

  full_file = argv[2];
  extract_file = argv[3];
  remainder_file = argv[4];

  /* Count how many items there are in total. */

  full = fopen(full_file,"r");
  if (full==NULL)
  { fprintf(stderr,"Can't open %s\n",full_file);
    exit(-1);
  }

  total = 0;
  ch = getc(full);
  while (ch!=EOF)
  { total += 1;
    do { ch = getc(full); } while (ch!='\n' && ch!=EOF);
    ch = getc(full);
  }
  
  fclose(full);

  /* Extract items. */

  full = fopen(full_file,"r");
  if (full==NULL)
  { fprintf(stderr,"Can't open %s\n",full_file);
    exit(-1);
  }

  extract = fopen(extract_file,"w");
  if (extract==NULL)
  { fprintf(stderr,"Can't create %s\n",extract_file);
    exit(-1);
  }

  remainder = fopen(remainder_file,"w");
  if (remainder==NULL)
  { fprintf(stderr,"Can't create %s\n",remainder_file);
    exit(-1);
  }

  rand_seed(seed);

  left = total;
  need = (int) (total*frac + 0.5);
  chosen = 0;
  ch = getc(full);

  while (ch!=EOF)
  { if (rand_uniform()<(double)(need-chosen)/left)
    { do { putc(ch,extract); ch = getc(full); } while (ch!='\n' && ch!=EOF);
      putc('\n',extract);
      chosen += 1;
    }
    else
    { do { putc(ch,remainder); ch = getc(full); } while (ch!='\n' && ch!=EOF);
      putc('\n',remainder);
    }
    left -= 1;
    ch = getc(full);
  }

  if (left!=0 || chosen!=need)
  { fprintf(stderr,
     "\nSOMETHING'S FISHY: total %d, need %d, chosen %d, left %d\n\n",
       total, need, chosen, left);
  }

  /* Print statistics. */

  fprintf(stderr, "%s: %d items out of %d, fraction %.3f aiming at %.3f\n",
    extract_file, chosen, total, (double)chosen/total, frac);

  exit(0);
}
