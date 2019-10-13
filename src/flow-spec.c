/* FLOW-SPEC.C - Specify flow model. */

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
  static flow_spec fspec, *flow = &fspec;  /* Static, so initialized to zeros */
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

    if ((flow = logg.data['F'])==0)
    { fprintf(stderr,"No flow specification in log file\n");
      exit(1);
    }

    printf("\n");

    if (flow->type=='t' || flow->type=='T')
    { 
      if (flow->type=='t')
      { printf("Flow model is \"test\"\n\n");
      }
      else
      { printf("Flow model is \"test-start-stop\"\n\n");
      }

      if (flow->lowU==flow->highU)
      { printf("Wind speed: %.3f\n",flow->lowU);
      }
      else
      { printf("Wind speed: %.3f:%.3f\n",flow->lowU,flow->highU);
      }

      printf("\nOther parameters:\n\n");
      printf("  ay: %.3f  by: %.5f  az: %.3f  bz: %.5f\n",
      flow->ay, flow->by, flow->az, flow->bz);
    }

    else /* type not recognized */
    { 
      printf("Unknown flow model: %c\n",flow->type);
      exit(1);
    }

    printf("\n");

    log_file_close(&logf);

    exit(0);
  }

  /* Otherwise, look at arguments. */

  if (argc<3) usage();

  if (strcmp(argv[2],"test")==0 || strcmp(argv[2],"test-start-stop")==0)
  {
    if (argc!=8) usage();

    flow->type = strcmp(argv[2],"test")==0 ? 't' : 'T';

    if (!parse_real_range(argv[3],&flow->lowU,&flow->highU,0)) usage();
    if (1/flow->highU==0) usage();
    if (flow->lowU<=0) 
    { fprintf(stderr,"Wind speed must be greater than zero\n");
      exit(1);
    }

    if (sscanf(argv[4],"%lf%c",&flow->ay,&junk)!=1) usage();
    if (sscanf(argv[5],"%lf%c",&flow->by,&junk)!=1) usage();
    if (sscanf(argv[6],"%lf%c",&flow->az,&junk)!=1) usage();
    if (sscanf(argv[7],"%lf%c",&flow->bz,&junk)!=1) usage();
  }

  else /* type not recognized */
  {
    fprintf(stderr,"Unknown type of flow model: %s\n",argv[2]);
    exit(1);
  }

  if (strchr(steady_state,flow->type))
  { if (src->max_start!=0 || src->max_stop!=1e30 || src->max_duration!=1e30)
    { fprintf(stderr,
        "Start/stop/duration times can't be used with steady-state models\n");
      exit(1);
    }
  }

  /* Append flow specification record to log file. */

  log_file_open(&logf,1);

  logf.header.type = 'F';
  logf.header.index = -1;
  logf.header.size = sizeof *flow;
  log_file_append(&logf,flow);

  log_file_close(&logf);

  exit(0);

}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
    "Usage: flow-spec log-file type { parameter }\n");
  fprintf(stderr,
    "   or: flow-spec log-file (to display stored specifications)\n");
  fprintf(stderr,
    "Parameters for \"test\" and \"test-start-stop\" models:  lowU[:highU] ay by az bz\n");

  exit(1);
}
