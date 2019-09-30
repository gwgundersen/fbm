/* LOG-LAST.C - Program to display index of last record in log files. */

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


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  log_file logf;
  int i;

  if (argc<2)
  { fprintf(stderr,"Usage: log-last log-file ...\n");
    exit(1);
  }

  for (i = 1; i<argc; i++)
  {
    logf.file_name = argv[i];
    log_file_open(&logf,0);
  
    log_file_last(&logf);
  
    if (logf.at_end)
    { printf("%s: Log file is empty\n",argv[i]);
    }
    else
    { printf("%s: Index of last record is %d\n",argv[i],logf.header.index);
    }

    log_file_close(&logf);
  }
  
  exit(0);
}
