/* MIX-SPEC.C - Program to specify a mixture model. */

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
#include "mix.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  static mix_spec spec, *mx = &spec;  /* Static so that unused fields are set */
                                      /*   to zero (for future compatibility) */
  log_file logf;
  log_gobbled logg;

  char ps[100];
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
    mix_record_sizes(&logg);

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }
  
    /* Display specifications. */
  
    printf("\n");
  
    if ((mx = logg.data['P'])==0)
    { printf ("No specification of mixture model priors found\n\n");
      exit(0);
    }

    printf ("Number of inputs:      %d\n", mx->N_inputs);  
    printf ("Number of targets:     %d\n", mx->N_targets);

    if (mx->N_components==0)
    { printf ("\nNumber of components:  Infinite\n");
    }
    else
    { printf ("\nNumber of components:  %d\n", mx->N_components);
    }

    printf("\n");

    printf ("Dirichlet concentration parameter: %s\n",
            prior_show(ps,mx->con_prior));
  
    printf ("Prior for SD hyperparameters:      %s\n",
            prior_show(ps,mx->SD_prior));
  
    printf ("Prior for mean component offsets:  %s\n",
            prior_show(ps,mx->mean_prior));
  
    printf("\n");
  
    log_file_close(&logf);
  
    exit(0);
  }

  /* Otherwise, figure out model form and priors from program arguments. */

  ap = argv+2;

  if (*ap==0 || (mx->N_inputs = atoi(*ap))<=0 && strcmp(*ap,"0")!=0) usage();
  ap += 1;

  if (*ap==0 || (mx->N_targets = atoi(*ap))<=0 && strcmp(*ap,"0")!=0) usage();
  ap += 1;

  mx->N_components = 0;

  if (*ap!=0 && strcmp(*ap,"/")!=0
       && (mx->N_components = atoi(*ap++))<=0) usage();

  if (*ap==0 || strcmp(*ap++,"/")!=0) usage();

  if (*ap==0 || !prior_parse(&mx->con_prior,*ap++))  usage();
  if (*ap==0 || !prior_parse(&mx->SD_prior,*ap++))   usage();

  if (*ap==0)
  { mx->mean_prior.width = 0;
    mx->mean_prior.scale = 0;
    mx->mean_prior.alpha[0] = 0;
    mx->mean_prior.alpha[1] = 0;
    mx->mean_prior.alpha[2] = 0;
  }
  else
  { if (!prior_parse(&mx->mean_prior,*ap++)) usage();
  }

  if (*ap!=0) usage();

  /* Check for illegal specifications. */

  if (mx->N_inputs!=0)
  { fprintf(stderr,"The number of inputs must be zero at present\n");
    exit(1);
  }

  if (mx->N_targets>Max_targets)
  { fprintf(stderr,"Too many targets (max %d)\n",Max_targets);
    exit(1);
  }

  if (mx->con_prior.alpha[0]!=0 
   || mx->con_prior.alpha[1]!=0 
   || mx->con_prior.alpha[2]!=0) 
  { fprintf(stderr,"Illegal concentration specification\n");
    exit(1);
  }

  if (mx->SD_prior.alpha[2]!=0)
  { fprintf(stderr,"Illegal prior for standard deviations\n");
    exit(1);
  }

  if (mx->mean_prior.scale
   || mx->mean_prior.alpha[0]!=0 
   || mx->mean_prior.alpha[1]!=0 
   || mx->mean_prior.alpha[2]!=0)
  { fprintf(stderr,"Illegal prior for means\n");
    exit(1);
  }

  if (mx->N_components==0 && !mx->con_prior.scale)
  { fprintf (stderr,
      "The prior for the concentration must be scaled when the number\n");
    fprintf (stderr,
      "  of components is infinite (ie, con-prior must start with \"x\")\n");
    exit(1);
  }

  /* Create log file and write records. */

  log_file_create(&logf);

  logf.header.type = 'P';
  logf.header.index = -1;
  logf.header.size = sizeof *mx;
  log_file_append(&logf,mx);

  log_file_close(&logf);

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
   "Usage: mix-spec log-file N-inputs N-targets [ N-components ]\n");
  fprintf(stderr,
   "                / concentration SD-prior [ mean-prior ]\n");
  fprintf(stderr,
   "       (Note: N-inputs must be 0 at present)\n");
  fprintf(stderr,
   "   or: mix-spec log-file (to display stored specifications)\n");
  exit(1);
}

