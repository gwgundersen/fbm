/* SRC-SPEC.C - Specify priors for source location model. */

/* Copyright (c) 2007 by Radford M. Neal 
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
#include "src.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
) 
{
  static src_spec sspec, *src = &sspec;  /* Static, so initialized to zeros */

  log_file logf;
  log_gobbled logg;
  char junk;
  int i, j;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing specifications. */

  if (argc==2)
  {
    /* Open log file and gobble up initial records. */
  
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    src_record_sizes(&logg);

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }

    /* Display specification. */  

    if ((src = logg.data['S'])==0)
    { fprintf(stderr,"No source specification in log file\n");
      exit(1);
    }

    printf("\n");

    if (src->highN==src->lowN)
    { printf("Number of sources: %d\n",src->lowN);
    }
    else
    { printf("Number of sources: %d:%d\n",src->lowN,src->highN);
    }

    if (src->highQ==src->lowQ && src->powQ==1)
    { printf("Source intensity:  %.2f\n",src->lowQ);
    }
    else if (src->powQ==1)
    { printf("Source intensity:  %.2f:%.2f\n",src->lowQ,src->highQ);
    }
    else
    { printf("Source intensity:  %.2f:%.2f:%.2f\n",
              src->lowQ,src->highQ,src->powQ);
    }

    if (src->max_start>0 || src->max_stop<1e30)
    { printf("\n");
    }
    if (src->max_start>0)
    { printf("Maximum start time: %.3f\n",src->max_start);
    }
    if (src->max_stop<1e30)
    { printf("Maximum stop time:  %.3f\n",src->max_stop);
    }
    if (src->max_duration<1e30)
    { printf("Maximum duration:   %.3f\n",src->max_duration);
    }

    printf("\nPrior ranges for coordinates of source locations:\n\n");
    for (j = 0; j<3 && (src->low[j]!=0 || src->high[j]!=0); j++)
    { printf("  %c: %.2f:%.2f","xyz"[j],src->low[j],src->high[j]);
    }
    printf("\n");
    
    printf("\n");

    log_file_close(&logf);

    exit(0);
  }

  /* Otherwise, look at arguments. */

  if (argc<5) usage();

  if (!parse_non_neg_range(argv[2],&src->lowN,&src->highN)) usage();
  if (!parse_real_range(argv[3],&src->lowQ,&src->highQ,&src->powQ)) usage();
  if (1/src->highQ==0) usage();
  if (src->lowQ<0) 
  { fprintf(stderr,"Source intensities must be non-negative\n");
    exit(1);
  }
  if (src->powQ<=0) 
  { fprintf(stderr,"Power for intensity transformation must be positive\n");
    exit(1);
  }

  src->max_start = 0;
  src->max_stop = 1e30;
  src->max_duration = 1e30;
  i = 4;

  if (argv[i]!=0 && strcmp(argv[i],"/")!=0 && *argv[i]!='+')
  { if (sscanf(argv[i],"%lf%c",&src->max_start,&junk)!=1) usage();
    if (src->max_start<0) 
    { fprintf(stderr,"Start times must be non-negative\n");
      exit(1);
    }
    i += 1;
  }
  if (argv[i]!=0 && strcmp(argv[i],"/")!=0 && *argv[i]!='+')
  { if (sscanf(argv[i],"%lf%c",&src->max_stop,&junk)!=1) usage();
    if (src->max_stop<src->max_start) 
    { fprintf(stderr,
         "Maximum stop time must not be less than maximum start time\n");
      exit(1);
    }
    i += 1;
  }
  if (argv[i]!=0 && *argv[i]=='+')
  { if (sscanf(argv[i],"%lf%c",&src->max_duration,&junk)!=1) usage();
    i += 1;
  }

  if (argv[i]==0 || strcmp(argv[i],"/")!=0) usage();
  i += 1;

  if (argv[i]==0) usage();

  for (j = 0; j<3 && argv[i+j]!=0; j++)
  { if (!parse_real_range(argv[i+j],&src->low[j],&src->high[j],0)) usage();
    if (1/src->high[j]==0) usage();
  }

  if (argv[i+j]!=0) usage();

  for ( ; j<3; j++)
  { src->low[j] = src->high[j] = 0;
  }

  if (src->low[2]<0)
  { fprintf(stderr,"The z (height) coordinate must be non-negative\n");
    exit(1);
  }

  /* Create log file and write records. */

  log_file_create(&logf);

  logf.header.type = 'S';
  logf.header.index = -1;
  logf.header.size = sizeof *src;
  log_file_append(&logf,src);

  log_file_close(&logf);

  exit(0);

}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr, "Usage: src-spec log-file\n");
  fprintf(stderr, "         lowN[:highN] lowQ[:highQ[:power]] [ max-start [ max-stop ] ]\n");
  fprintf(stderr, "           / lowx:highx [ lowy:highy [ lowz:highz ] ]\n");
  fprintf(stderr,
                  "   or: src-spec log-file (to display stored specifications)\n");

  exit(1);
}
