/* RAND-SEED.C - Program for specifying random number seed. */

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

#include "rand.h"
#include "log.h"


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{ 
  log_file logf;
  log_gobbled logg;

  rand_state *rs;

  int seed;

  if (argc<2 || argc>3 
   || argc==3 && (seed = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0)
  {
    fprintf(stderr,"Usage: rand-seed log-file seed\n");
    fprintf(stderr,"   or: rand-seed log-file (to display stored seed)\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* See if we are to display existing seed. */

  if (argc==2)
  {
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    logg.req_size['r'] = sizeof *rs;

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }

    if (logg.data['r']==0)
    { printf("\nNo random number seed stored in log file with index -1\n\n");
    }
    else
    { rs = logg.data['r'];
      printf("\nRandom number seed: %d\n\n",rs->seed);
    }
   
    exit(0);
  }

  /* Otherwise, append state structure initialized using given seed. */
  
  rand_seed(seed);

  rs = rand_get_state();

  log_file_open (&logf, 1);

  log_file_last (&logf);
  if (logf.at_end)
  { fprintf(stderr,"Log file is empty\n");
    exit(1);
  }

  logf.header.type = 'r';
  logf.header.size = sizeof *rs;
  log_file_append (&logf, rs);

  log_file_close (&logf);

  exit(0);

}
