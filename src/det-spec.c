/* DET-SPEC.C - Specify detector noise model. */

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

int main
( int argc,
  char **argv
) 
{
  static det_spec dspec, *det = &dspec;  /* Static, so initialized to zeros */
  src_spec *src;

  log_file logf;
  log_gobbled logg;
  char junk;
  int i, j;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* Open log file and gobble up initial records. */
  
  log_file_open(&logf,0);

  log_gobble_init(&logg,0);
  src_record_sizes(&logg);

  if (!logf.at_end && logf.header.index==-1)
  { log_gobble(&logf,&logg);
  }

  if ((src = logg.data['S'])==0)
  { fprintf(stderr,"No source specification in log file\n");
    exit(1);
  }

  /* See if we are to display existing specifications. */

  if (argc==2)
  {
    /* Display specification. */  

    if ((det = logg.data['T'])==0)
    { fprintf(stderr,"No detector specification in log file\n");
      exit(1);
    }

    printf("\n");

    if (det->log_low_width==det->log_high_width)
    { printf("Width of noise distribution: %.4f\n",exp(det->log_low_width));
    }
    else
    { printf("Width of noise distribution: %.4f:%.4f\n",
              exp(det->log_low_width), exp(det->log_high_width));
    }

    if (det->inv_low_df==0)
    { printf("\nNoise has Gaussian distribution\n");
    }
    else 
    { printf("\nNoise has t distribution\n\n");
      if (det->inv_low_df==det->inv_high_df)
      { printf("Degrees of freedom for noise: %.2f\n",1/det->inv_low_df);
      }
      else if (det->inv_high_df==0)
      { printf("Degrees of freedom for noise: %.2f:\n",1/det->inv_low_df);
      }
      else
      { printf("Degrees of freedom for noise: %.2f:%.2f\n",
                1/det->inv_low_df, 1/det->inv_high_df);
      }
    }
    
    printf("\n");

    log_file_close(&logf);

    exit(0);
  }

  /* Otherwise, look at arguments. */

  if (argc<3) usage();

  if (!parse_real_range(argv[2],&det->log_low_width,&det->log_high_width,0)
   || det->log_low_width<=0 || 1/det->log_high_width==0)
  { usage();
  }
  det->log_low_width = log(det->log_low_width);
  det->log_high_width = log(det->log_high_width);

  if (argc>=4)
  { if (!parse_real_range(argv[3],&det->inv_low_df,&det->inv_high_df,0)) 
    { usage();
    }
    if (det->inv_low_df<=0) usage();
    det->inv_low_df = 1/det->inv_low_df;
    det->inv_high_df = 1/det->inv_high_df;
  }

  if (argc>4) usage();

  /* Append detector specification record to log file. */

  log_file_open(&logf,1);

  logf.header.type = 'T';
  logf.header.index = -1;
  logf.header.size = sizeof *det;
  log_file_append(&logf,det);

  log_file_close(&logf);

  exit(0);

}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
    "Usage: det-spec log-file low-width[:high-width] [ low-df[:[high-df]] ]\n");
  fprintf(stderr,
    "   or: det-spec log-file (to display stored specifications)\n");

  exit(1);
}
