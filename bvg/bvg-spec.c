/* BVG-SPEC.C - Specify bivariate Gaussian. */

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
#include "mc.h"
#include "bvg.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
) 
{
  bvg_spec bspec, *bs = &bspec;

  log_file logf;
  log_gobbled logg;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing specifications. */

  if (argc==2)
  {
    /* Open log file and gobble up initial records. */
  
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    logg.req_size['B'] = sizeof *bs;

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }

    /* Display specification. */  

    if ((bs = logg.data['B'])==0)
    { fprintf(stderr,"No bivariate Gaussian specification in log file\n");
      exit(1);
    }

    printf("\n");
    printf("Standard deviations: %.6f %.6f\n",bs->std1,bs->std2);
    printf("Correlation:         %.6f\n",bs->corr);
    if (bs->rep>1) printf("Replications:        %d\n",bs->rep);
    printf("\n");

    log_file_close(&logf);

    exit(0);
  }

  /* Otherwise, look at arguments. */

  if (argc!=5 && argc!=6) usage();

  bs->rep = 1;

  if ((bs->std1 = atof(argv[2]))<=0
   || (bs->std2 = atof(argv[3]))<=0
   || (bs->corr = atof(argv[4]))==0 && strcmp(argv[4],"0")!=0
   || bs->corr<=-1.0 || bs->corr>=1.0
   || argc==6 && (bs->rep = atoi(argv[5]))<=0)
  { usage();
  }

  /* Create log file and write records. */

  log_file_create(&logf);

  logf.header.type = 'B';
  logf.header.index = -1;
  logf.header.size = sizeof *bs;
  log_file_append(&logf,bs);

  log_file_close(&logf);

  exit(0);

}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
    "Usage: bvg-spec log-file std1 std2 corr [ replication ]\n");
  fprintf(stderr,
    "   or: bvg-spec log-file (to display stored specifications)\n");

  exit(1);
}
