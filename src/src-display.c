/* SRC-DISPLAY.C - Print parameters of source location model. */

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


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  src_spec *src;
  det_spec *det;
  flow_spec *flow;
  src_params *params;

  log_file logf;
  log_gobbled logg;

  int index;
  char junk;
  int i, j;

  /* Look at arguments. */

  index = -1;

  if (argc!=2 && argc!=3 || argc>2 && sscanf(argv[2],"%d%c",&index,&junk)!=1)
  { fprintf (stderr, "Usage: src-display log-file [ index ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  src_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  src = logg.data['S'];
  det = logg.data['T'];
  flow = logg.data['F'];

  src_check_specs_present (src, det, flow);

  /* Read the desired state from the log file. */

  if (index<0)
  { 
    log_gobble_last(&logf,&logg);

    if (logg.last_index<0)
    { fprintf(stderr,"No source location parameters in log file\n");
      exit(1);
    }

    index = logg.last_index;
  }
  else
  {
    while (!logf.at_end && logf.header.index!=index)
    { log_file_forward(&logf);
    }

    if (logf.at_end)
    { fprintf(stderr,"That index does not appear in the log file\n");
      exit(1);
    }

    log_gobble(&logf,&logg);
  }

  if (logg.index['q']!=index)
  { fprintf(stderr,"No source location parameters stored with that index\n");
    exit(1);
  }

  params = logg.data['q'];

  printf ("\nPARAMETERS OF SOURCE LOCATION MODEL AT ITERATION %d\n",index);

  /* Print parameters of noise model. */

  if (det->inv_low_df==0)
  { printf ("\nNoise standard deviation: %.4f\n", 
             exp(params->log_width));
  }
  else 
  { printf ("\nNoise width: %.4f, df: %.3f\n", 
             exp(params->log_width), 1/params->inv_df);
  }

  /* Print parameters of flow model (if any). */

  if (flow->type=='t')
  { printf ("\nWind speed: %.4f\n", params->U);
  }

  /* Print the parameters of each source. */

  if (strchr(steady_state,flow->type))
  { printf (
      "\n        Q       x       y       z\n\n");
  }
  else if (src->max_stop==1e30 && src->max_duration==1e30)
  { printf (
      "\n        Q       x       y       z     start\n\n");
  }
  else
  { printf (
      "\n        Q       x       y       z     start    stop  duration\n\n");
  }

  for (i = 0; i<(int)params->N0; i++)
  { 
    if (strchr(steady_state,flow->type))
    { printf("%3d %7.2f %7.2f %7.2f %7.2f\n", i+1,
        pow(params->src[i].Q,1/src->powQ),
        params->src[i].coord[0], 
        params->src[i].coord[1], 
        params->src[i].coord[2]);
    }
    else if (src->max_stop==1e30 && src->max_duration==1e30)
    { printf("%3d %7.2f %7.2f %7.2f %7.2f %7.2f\n", i+1,
        pow(params->src[i].Q,1/src->powQ),
        params->src[i].coord[0], 
        params->src[i].coord[1], 
        params->src[i].coord[2],
        params->src[i].start);
    }
    else
    { printf("%3d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n", i+1,
        pow(params->src[i].Q,1/src->powQ),
        params->src[i].coord[0], 
        params->src[i].coord[1], 
        params->src[i].coord[2],
        params->src[i].start, 
        params->src[i].stop,
        params->src[i].stop - params->src[i].start);
    }
  }
  
  exit(0);
}
