/* SRC-INITIAL.C - Program to set initial state of source location model. */

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
#include "rand.h"

static void usage();


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  log_file logf;
  log_gobbled logg;

  src_spec *src;
  det_spec *det;
  flow_spec *flow;

  src_params *params; 

  char junk;
  double v;
  int a, i, j;
 
  /* Get name of log file. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  src_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  src = logg.data['S'];
  det = logg.data['T'];
  flow = logg.data['F'];

  src_check_specs_present (src, det, flow);

  /* Allocate space for parameters.  Allow enough space for the specified 
     maximum number of sources. */

  params = 
    chk_alloc (src_params_length(src->highN)/sizeof(double), sizeof(double));

  /* Set parameters to defaults to start. */

  src_default_parameters (src,det,flow,params);

  /* Set flow parameters according to command arguments. */

  a = 2;

  if (flow->type=='t')
  { if (argv[a]!=0 && strcmp(argv[a],"/")!=0)
    { if (sscanf(argv[a],"%lf%c",&params->U,&junk)!=1) usage();
      if (params->U<flow->lowU || params->U>flow->highU)
      { fprintf(stderr,"Wind speed out of range\n");
        exit(1);
      }
      a += 1;
    }
  }

  if (argv[a]==0 || strcmp(argv[a],"/")!=0) usage();

  a += 1;

  /* Set detector parameters according to command arguments. */

  if (argv[a]!=0 && strcmp(argv[a],"/")!=0)
  { 
    if (sscanf(argv[a],"%lf%c",&v,&junk)!=1 || v<0) usage();
    params->log_width = log(v);
    if (params->log_width<det->log_low_width 
     || params->log_width>det->log_high_width)
    { fprintf(stderr,"Detector noise width out of range\n");
      exit(1);
    }
 
    a += 1;

    if (argv[a]!=0 && strcmp(argv[a],"/")!=0) 
    {
      if (sscanf(argv[a],"%lf%c",&v,&junk)!=1 || v<0) usage();
      params->inv_df = 1/v;
      if (params->inv_df>det->inv_low_df || params->inv_df<det->inv_high_df)
      { fprintf(stderr,"Detector noise degrees of freedom out of range\n");
        exit(1);
      }
      a += 1;
    }
  }

  /* Set source parameters according to command arguments. */

  i = 0;

  while (argv[a]!=0)
  { if (strcmp(argv[a],"/")!=0) usage();
    a += 1;

    if (i==src->highN)
    { fprintf(stderr,"Initial values for too many sources specified\n");
      exit(1);
    }

    if (argv[a]==0 || sscanf(argv[a],"%lf%c",&v,&junk)!=1 || v<0) usage();
    if (v<src->lowQ || v>src->highQ)
    { fprintf(stderr,"Source intensity out of range\n");
      exit(1);
    }
    params->src[i].Q = pow(v,src->powQ);
    a += 1;

    for (j = 0; j<3; j++)
    {
      if (argv[a]==0 || sscanf(argv[a],"%lf%c",&v,&junk)!=1) usage();
      if (v<src->low[j] || v>src->high[j])
      { fprintf(stderr,"Source location out of range\n");
        exit(1);
      }
      params->src[i].coord[j] = v;
      a += 1;
    }

    if (argv[a]!=0 && strcmp(argv[a],"/")!=0)
    { 
      if (sscanf(argv[a],"%lf%c",&v,&junk)!=1 || v<0) usage();
      if (v<0 || v>src->max_start)
      { fprintf(stderr,"Start time out of range\n");
        exit(1);
      }
      params->src[i].start = v;
      a += 1;

      if (argv[a]!=0 && strcmp(argv[a],"/")!=0)
      { 
        if (src->max_stop>=1e30 && src->max_duration==1e30)
        { fprintf(stderr,"A stop time should not be specified\n");
          exit(1);
        }
        if (sscanf(argv[a],"%lf%c",&v,&junk)!=1 || v<0) usage();
        if (v<params->src[i].start || v>src->max_stop  
         || v-params->src[i].start>src->max_duration)
        { fprintf(stderr,"Stop time out of range\n");
          exit(1);
        }
        params->src[i].stop = v;
        a += 1;
      }
    }
    
    i += 1;
  }

  params->N0 = i>src->lowN ? i : src->lowN;

  /* Write out parameter record. */

  logf.header.type = 'q';
  logf.header.index = 0;
  logf.header.size = src_params_length(src->highN);
  log_file_append (&logf, params);

  log_file_close(&logf);

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr, 
"Usage: src-initial log-file { flow-parameter } / [ noise-width [ noise-df ] ]\n");
  fprintf(stderr,
"                            { / Q x y z [ start [ stop ] ] }\n");

  exit(1);
}
