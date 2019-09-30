/* GP-EIGEN.C - Program to find eigenvalues/vectors of Gaussian process. */

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

  double *cov, *eigen, *evec, *tvec;
  double *inputs;
  int N_cases;

  int vec, val;

  char **ap;
  double v;
  int i, j, k;

  /* Look at arguments. */

  vec = 0; val = 0;

  while (argc>1)
  { if (strcmp(argv[1],"-values")==0)
    { val = 1;
    }
    else if (strcmp(argv[1],"-vectors")==0)
    { vec = 1;
    }
    else
    { break;
    }
    argc -= 1;
    argv += 1;
  }

  if (vec==0 && val==0) 
  { val = 1;
  }
  
  if (argc<3) usage();

  logf.file_name = argv[1];
  
  if ((index = atoi(argv[2]))<0) usage();

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

  if (extra == 0 && (m!=0 && m->noise.alpha[2]!=0))
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
  { strcpy(data_spec->test_inputs,input_file);
    gp_data_read (0, 1, gp, m, 0);
    N_cases = N_test;
    inputs = test_inputs;
  }
  else
  { gp_data_read (1, 0, gp, m, 0);
    N_cases = N_train;
    inputs = train_inputs;
  }

  /* Allocate space for covariance matix and eigenvalues/vectors. */

  cov = chk_alloc (N_cases*N_cases, sizeof (double));
  eigen = chk_alloc (N_cases, sizeof (double));
  if (vec)   
  { evec = chk_alloc (N_cases*N_cases, sizeof (double));
    tvec = chk_alloc (N_cases, sizeof (double)); 
  }

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

  /* Compute covariance for latent values at training points.*/

  gp_cov(gp, h, inputs, N_cases, inputs, N_cases, cov, 0, 0);

  if (gp->has_jitter)
  { for (i = 0; i<N_cases; i++)
    { cov[i*N_cases+i] += exp(2 * *h->jitter);
    }
  }

  /* Now add on extra, as requested. */

  if (extra>0)
  { for (i = 0; i<N_cases; i++) 
    { cov[i*N_cases+i] += extra;
    }
  }
  else if (extra==0)
  { if (m==0 || m->type!='R' || m->noise.alpha[2]!=0) abort();
    for (i = 0; i<N_cases; i++) 
    { cov[i*N_cases+i] += exp(2 * *h->noise[0]);
    }
  }

  /* Find eigenvalues of covariance matrix. */

  jacobi (cov, vec ? evec : 0, 1e-10, N_cases);

  /* Sort eigenvalues. */

  for (i = 0; i<N_cases; i++)
  { v = cov[i*N_cases+i];
    if (vec) 
    { for (k = 0; k<N_cases; k++) tvec[k] = evec[i*N_cases+k];
    }
    for (j = i; j>0 && eigen[j-1]<v; j--) 
    { eigen[j] = eigen[j-1];
      if (vec) 
      { for (k = 0; k<N_cases; k++) evec[j*N_cases+k] = evec[(j-1)*N_cases+k];
      }
    }
    eigen[j] = v;
    if (vec)
    { for (k = 0; k<N_cases; k++) evec[j*N_cases+k] = tvec[k];
    }
  }

  /* And finally, print eigenvalues and/or eigenvectors. */

  for (i = 0; i<N_cases; i++)
  { if (val) 
    { printf("%.8e",eigen[i]);
    }
    if (vec)
    { for (k = 0; k<N_cases; k++)
      { printf(" %+.8e",evec[i*N_cases+k]);
      }
    }
    printf("\n");
  }

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
    "Usage: gp-eigen [ -values ] [ -vectors ] log-file index [ extra ]\n");
  fprintf(stderr,
    "                [ / test-file ]\n");
  exit(1);
}
