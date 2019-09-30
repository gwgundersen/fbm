/* LOG-LAST.C - Program to display index of last record in log file. */

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

  if (argc!=2)
  { fprintf(stderr,"Usage: log-last log-file\n");
    exit(1);
  }

  logf.file_name = argv[1];
  log_file_open(&logf,0);

  log_file_last(&logf);

  if (logf.at_end)
  { printf("Log file is empty\n");
  }
  else
  { printf("Index of last record: %d\n",logf.header.index);
  }
  
  exit(0);
}
