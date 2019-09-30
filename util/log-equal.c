/* LOG-EQUAL.C - Check if records at given indexes in two log files match. */

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

#include "log.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  char *ignore;
  log_file logf1, logf2;
  int ix1, ix2;
  log_gobbled logg1, logg2;
  char junk;
  int i;

  ignore = "";
  if (argc>1 && *argv[1]=='-')
  { ignore = argv[1]+1;
    argc -= 1;
    argv += 1;
  }

  if (argc!=5) usage();

  if (sscanf(argv[2],"%d%c",&ix1,&junk)!=1
   || sscanf(argv[4],"%d%c",&ix2,&junk)!=1)
  { usage();
  }

  logf1.file_name = argv[1];
  log_file_open(&logf1,0);

  logf2.file_name = argv[3];
  log_file_open(&logf2,0);

  while (!logf1.at_end && logf1.header.index<ix1)
  { log_file_forward(&logf1);
  }

  if (logf1.at_end || logf1.header.index!=ix1)
  { fprintf(stderr,"No records at index %d of file %s\n",ix1,logf1.file_name);
    exit(-1);
  }

  while (!logf2.at_end && logf2.header.index<ix2)
  { log_file_forward(&logf2);
  }

  if (logf2.at_end || logf2.header.index!=ix2)
  { fprintf(stderr,"No records at index %d of file %s\n",ix2,logf2.file_name);
    exit(-1);
  }

  log_gobble_init(&logg1,0);
  log_gobble(&logf1,&logg1);

  log_gobble_init(&logg2,0);
  log_gobble(&logf2,&logg2);

  for (i = 1; i<128; i++)
  { if (!strchr(ignore,i))
    { if (logg1.index[i]==ix1 || logg2.index[i]==ix2)
      { if (logg1.index[i]!=ix1 || logg2.index[i]!=ix2
         || logg1.actual_size[i]!=logg2.actual_size[i]
         || memcmp(logg1.data[i],logg2.data[i],logg1.actual_size[i])!=0)
        { exit(1);
        }
      }
    }
  }

  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage(void)
{ 
  fprintf(stderr,
"Usage: log-equal [ -ignoredtype... ] log-file-1 index-1 log-file-2 index-2\n");
  exit(1);
}
