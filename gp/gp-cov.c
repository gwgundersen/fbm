/* GP-COV.C - Program to print covariance matrix for Gaussian process. */

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
#include "data.h"
#include "prior.h"
#include "model.h"
#include "matrix.h"
#include "rand.h"
#include "gp.h"
#include "gp-data.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  gp_spec *gp;
  model_specification *m;
  gp_hypers hypers, *h = &hypers;

  log_file logf;
  log_gobbled logg;

  double extra;
  char *input_file;
  int index;

  double *cov;

  char **ap;
  int i, j;

  /* Look at arguments. */
  
  if (argc<3) usage();

  logf.file_name = argv[1];
  
  if ((index = atoi(argv[2]))<0 || index==0 && strcmp(argv[2],"0")!=0) usage();

  extra = -1;
 
  ap = argv+3;

  if (*ap!=0 && strcmp(*ap,"/")!=0)
  { if (strcmp(*ap,"+noise")==0)
    { extra = 0;
    }
    else 
    { if ((extra = atof(*ap))<=0) usage();
    }
    ap += 1;
  }

  input_file = 0;
  
  if (*ap!=0)
  { if (strcmp(*ap++,"/")!=0) usage();
    input_file = *ap++;
  }

  if (*ap!=0) usage();

  /* Open log file and read Gaussian process specifications. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  gp_record_sizes(&logg);

  if (!logf.at_end && logf.header.index==-1)
  { log_gobble(&logf,&logg);
  }

  gp = logg.data['P'];
  m = logg.data['M'];

  if (gp==0)
  { fprintf(stderr,"No specification for Gaussian process in log file\n");
    exit(1);
  }

  if (extra==0 && (m==0 || m->type!='R'))
  { fprintf(stderr,"Can't add noise if the model is not for real data\n");
    exit(1);
  }

  if (extra==0 && (m!=0 && m->noise.alpha[2]!=0))
  { fprintf(stderr,"Can't add noise when it varies from case to case\n");
    exit(1);
  }

  /* Read data specification. */

  data_spec = logg.data['D'];

  if (data_spec==0)
  { fprintf(stderr,"No data specification\n");
    exit(1);
  }

  if (input_file) 
  { strcpy(data_spec->train_inputs,input_file);
  }

  data_spec->train_targets[0] = 0; /* Don't need the target values */

  gp_data_read (1, 0, gp, m, 0);

  /* Allocate space for covariance matix. */

  cov = chk_alloc (N_train*N_train, sizeof (double));

  /* Skip to desired index. */

  while (!logf.at_end && logf.header.index<index)
  { log_file_forward(&logf);
  }

  if (logf.at_end || logf.header.index>index)
  { fprintf(stderr,"No records stored with specified index\n");
    exit(1);
  }
  
  /* Read the desired network from the log file. */

  h->total_hypers = gp_hyper_count(gp,m);
  logg.req_size['S'] = h->total_hypers * sizeof(double);

  log_gobble(&logf,&logg);
       
  if (logg.data['S']==0 || logg.index['S']!=logg.last_index)
  { fprintf(stderr, "Record missing at index %d\n",logg.last_index);
    exit(1);
  }
  
  h->hyper_block = logg.data['S'];
  gp_hyper_pointers (h, gp, m);

  /* Compute covariance for latent values or targets at training points.*/

  gp_train_cov (gp, m, h, 0, 0, extra!=0 ? cov : 0, extra==0 ? cov : 0, 0);

  /* Now add on extra, if requested. */

  if (extra>0)
  { for (i = 0; i<N_train; i++) 
    { cov[i*N_train+i] += extra;
    }
  }

  /* And finally print covariances. */

  for (i = 0; i<N_train; i++)
  { for (j = 0; j<N_train; j++)
    { printf(" %.18e",cov[i*N_train+j]);
    }
    printf("\n");
  }

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,"Usage: gp-cov log-file index [ extra ] [ / train-inputs ]\n");
  exit(1);
}
