/* LOG-RECORDS.C - Program to list records present in log file. */

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

void main
( int argc,
  char **argv
)
{
  log_file logf;
  int index;

  if (argc!=2)
  { fprintf(stderr,"Usage: log-records log-file\n");
    exit(1);
  }

  logf.file_name = argv[1];
  log_file_open(&logf,0);

  while (!logf.at_end)
  {
    index = logf.header.index;
    printf("%5d",index);

    while (!logf.at_end && logf.header.index==index)
    { printf(" %c:%d",logf.header.type,logf.header.size);
      log_file_forward(&logf);
    }

    printf("\n");
  }
  
  exit(0);
}
