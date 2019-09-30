/* LOG-APPEND.C - Append records from one log file to another log file. */

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

#include "misc.h"
#include "log.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  log_file logf_in, logf_out;
  int low, high, modulus;
  int last, index, position;
  char junk;
  void *data;
  char *ignore;

  ignore = "";
  if (argc>1 && *argv[1]=='-')
  { ignore = argv[1]+1;
    argc -= 1;
    argv += 1;
  }

  position = -1;

  if (strcmp(argv[argc-1],"-")==0 
   || sscanf(argv[argc-1],"%d%c",&position,&junk)==1)
  { 
    if (strcmp(argv[argc-1],"-")==0)
    { position = -2;
    }

    argc -=1 ;
  }

  if (argc!=3 && argc!=4) usage();

  logf_in.file_name = argv[1];
  log_file_open(&logf_in,0);
  log_file_last(&logf_in);

  if (logf_in.at_end)
  { fprintf(stderr,"Input log file is empty\n");
    exit(1);
  }
  
  if (logf_in.header.index<0)
  { fprintf(stderr,"Log file has no records at positive indexes\n");
    exit(1);
  }

  if (argc==4)
  { parse_range(argv[2],&low,&high,&modulus);
    if (modulus<0)
    { fprintf(stderr,"Bad range specification: %s\n",argv[2]);
      exit(1);
    }
    if (high==-2) high = low;
  }
  else
  { low = high = logf_in.header.index;
    modulus = 1;
  }

  logf_out.file_name = argv[argc-1];
  log_file_open(&logf_out,1);
  log_file_last(&logf_out);

  index = position>=0 ? position-1
        : logf_out.at_end || logf_out.header.index<0 ? -1 
        : position==-2 ? logf_out.header.index-1 : logf_out.header.index; 

  last = -1;

  log_file_first(&logf_in);

  while (!logf_in.at_end && (high<0 || logf_in.header.index<=high))
  {
    data = 0;

    if (!strchr(ignore,logf_in.header.type) &&
        logf_in.header.index>=low && logf_in.header.index%modulus==0)
    { 
      data = malloc(logf_in.header.size);
      if (data==0)
      { fprintf(stderr,"Not enough memory for record of size %d - ignored\n",
          logf_in.header.size);
      }
    }

    if (data==0)
    { log_file_forward(&logf_in);
    }
    else
    {
      if (logf_in.header.index!=last)
      { last = logf_in.header.index;
        index += 1;
      }

      logf_out.header = logf_in.header;
      logf_out.header.index = index;

      log_file_read (&logf_in, data, logf_in.header.size);

      log_file_append (&logf_out, data);

      free(data);
    } 
  }

  log_file_close(&logf_in);
  log_file_close(&logf_out);
  
  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage(void)
{ 
  fprintf(stderr,
  "Usage: log-append [ -ignoredtype... ] logfile-in [ range ] logfile-out [ index ]\n");
  exit(1);
}
