/* DFT-SPEC.C - Program to specify a Dirichlet diffusion tree model. */

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
#include "prior.h"
#include "model.h"
#include "data.h"
#include "dft.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  static dft_spec spec, *dft = &spec; /* Static so that unused fields are set */
                                      /*   to zero (for future compatibility) */
  log_file logf;
  log_gobbled logg;

  char ps[100];
  int dt;
  char **ap;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing specifications. */

  if (argc==2)
  {
    /* Open log file and gobble up initial records. */
  
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    dft_record_sizes(&logg);

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }
  
    /* Display specifications. */
  
    printf("\n");
  
    if ((dft = logg.data['P'])==0)
    { printf ("No specification of diffusion tree model found\n\n");
      exit(0);
    }

    printf ("Number of inputs:  %d\n", dft->N_inputs);  
    printf ("Number of targets: %d\n", dft->N_targets);
    printf ("Number of trees:   %d\n", dft->N_trees);

    for (dt = 0; dt<dft->N_trees; dt++)
    { printf("\nTree %d:\n\n",dt+1);
      printf ("  Prior for standard deviations: %s\n",
               prior_show(ps,dft->tree[dt].d_SD));
      printf ("  Priors for divergence function parameters:");
      printf (" %s",  dft->tree[dt].c0.width==0 ? " -" :
                      prior_show(ps,dft->tree[dt].c0));
      printf (" %s",  dft->tree[dt].c1.width==0 ? " -" :
                      prior_show(ps,dft->tree[dt].c1));
      printf (" %s",  dft->tree[dt].c2.width==0 ? " -" :
                      prior_show(ps,dft->tree[dt].c2));
      printf("\n");
    }

    printf("\n");
  
    log_file_close(&logf);
  
    exit(0);
  }

  /* Otherwise, figure out model and priors from program arguments. */

  ap = argv+2;

  if (*ap==0 || (dft->N_inputs = atoi(*ap))<=0 && strcmp(*ap,"0")!=0) usage();
  ap += 1;

  if (*ap==0 || (dft->N_targets = atoi(*ap))<=0 && strcmp(*ap,"0")!=0) usage();
  ap += 1;

  dft->N_trees = 0;

  while (*ap!=0 && strcmp(*ap++,"/")==0) 
  {
    if (dft->N_trees==Max_trees)
    { fprintf(stderr,"Too many diffusion trees (max %d)\n",Max_trees);
      exit(1);
    }

    if (*ap==0 || !prior_parse(&dft->tree[dft->N_trees].d_SD,*ap++)) usage();

    if (*ap!=0 && strcmp(*ap,"/")!=0)
    { if (strcmp(*ap,"-")!=0 && strcmp(*ap,"0")!=0)
      { if (!prior_parse(&dft->tree[dft->N_trees].c0,*ap)) usage();
      }
      ap += 1;
    }

    if (*ap!=0 && strcmp(*ap,"/")!=0)
    { if (strcmp(*ap,"-")!=0 && strcmp(*ap,"0")!=0)
      { if (!prior_parse(&dft->tree[dft->N_trees].c1,*ap)) usage();
      }
      ap += 1;
    }

    if (*ap!=0 && strcmp(*ap,"/")!=0)
    { if (strcmp(*ap,"-")!=0 && strcmp(*ap,"0")!=0)
      { if (!prior_parse(&dft->tree[dft->N_trees].c2,*ap)) usage();
      }
      ap += 1;
    }

    dft->N_trees += 1;
  }

  if (dft->N_trees==0)
  { fprintf(stderr,"Model must have at least one diffusion tree\n");
    exit(1);
  }

  if (*ap!=0) usage();

  /* Check for illegal specifications. */

  if (dft->N_inputs!=0)
  { fprintf(stderr,"The number of inputs must be zero at present\n");
    exit(1);
  }

  if (dft->N_targets>Max_targets)
  { fprintf(stderr,"Too many targets (max %d)\n",Max_targets);
    exit(1);
  }

  for (dt = 0; dt<dft->N_trees; dt++)  
  { if (dft->tree[dt].d_SD.scale || dft->tree[dt].d_SD.alpha[2]!=0)
    { fprintf(stderr,"Illegal prior for standard deviations\n");
      exit(1);
    }
    if (dft->tree[dt].c0.scale || dft->tree[dt].c0.alpha[2]!=0 
                               || dft->tree[dt].c0.alpha[1]!=0)
    { fprintf(stderr,"Illegal prior for c0\n");
      exit(1);
    }
    if (dft->tree[dt].c1.scale || dft->tree[dt].c1.alpha[2]!=0 
                               || dft->tree[dt].c1.alpha[1]!=0)
    { fprintf(stderr,"Illegal prior for c1\n");
      exit(1);
    }
    if (dft->tree[dt].c2.scale || dft->tree[dt].c2.alpha[2]!=0 
                               || dft->tree[dt].c2.alpha[1]!=0)
    { fprintf(stderr,"Illegal prior for c2\n");
      exit(1);
    }
  }

  /* Create log file and write records. */

  log_file_create(&logf);

  logf.header.type = 'P';
  logf.header.index = -1;
  logf.header.size = sizeof *dft;
  log_file_append(&logf,dft);

  log_file_close(&logf);

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
"Usage: dft-spec log-file N-inputs N-targets { / SD-prior c0 [ c1 [ c2 ] ] }\n"
);
  fprintf(stderr,
"       (Note: N-inputs must be 0 at present)\n");
  fprintf(stderr,
"   or: dft-spec log-file (to display stored specifications)\n");
  exit(1);
}
