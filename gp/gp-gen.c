/* GP-GEN.C - Program to generate Gaussian process hyperparameters. */

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
#include "matrix.h"
#include "log.h"
#include "prior.h"
#include "data.h"
#include "model.h"
#include "gp.h"
#include "gp-data.h"
#include "rand.h"


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  gp_spec *gp;

  model_specification *m;
  model_survival *v;

  gp_hypers hypers, *h = &hypers;

  int gen_latent, gen_noise_var, max_index, index, fixarg, fix;
  double scale_value, relevance_value;
  double *latent_values, *latent_cov, *scratch, *noise_variances;

  log_file logf;
  log_gobbled logg;
 
  /* Look at program arguments. */

  max_index = 0;
  scale_value = 0;
  relevance_value = 0;
  gen_latent = 0;
  gen_noise_var = 0;

  while (argc>1)
  { 
    if (strcmp(argv[1],"-l")==0)
    { gen_latent = 1;
    }
    else if (strcmp(argv[1],"-n")==0)
    { gen_noise_var = 1;
    }
    else
    { break;
    }

    argc -= 1;
    argv += 1;
  }

  fixarg = argc>2 && strcmp(argv[2],"fix")==0 ? 2
         : argc>3 && strcmp(argv[3],"fix")==0 ? 3
         : argc;

  fix = fixarg!=argc;

  if (argc<2 || argc!=fixarg && argc!=fixarg+1 && argc!=fixarg+3
   || fixarg>2 && (max_index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0
   || fixarg+1<argc && (scale_value = atof(argv[fixarg+1]))<=0
   || fixarg+2<argc && (relevance_value = atof(argv[fixarg+2]))<=0)
  { fprintf(stderr,
"Usage: gp-gen [ -l ] [ -n ] log-file [ max-index ] [ \"fix\" [ sc-val rel-val ] ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  gp_record_sizes(&logg);
  logg.req_size['r'] = sizeof (rand_state);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  gp = logg.data['P'];
  m = logg.data['M'];
  v = logg.data['V'];

  gp_check_specs_present(gp,0,m,v);

  if (gen_noise_var)
  { if (m==0 || m->type!='R' || m->noise.alpha[2]==0)
    { gen_noise_var = 0;  /* Silently ingore the -n option */
    }
  }

  /* Read training inputs if generating latent values and/or noise variances. */

  if (gen_latent || gen_noise_var)
  {  
    if (m!=0 && m->type=='V')
    { fprintf(stderr,"Can't handle survival models in gp-gen\n");
    }

    data_spec = logg.data['D'];

    if (data_spec==0)
    { fprintf (stderr,
"Can't generate latent values or noise variances if there's no data specification\n");
      exit(1);
    }

    if (logg.actual_size['D'] !=
          data_spec_size(data_spec->N_inputs,data_spec->N_targets))
    { fprintf(stderr,"Data specification record is the wrong size!\n");
      exit(1);
    }

    data_spec->train_targets[0] = 0;  /* Don't need targets for training cases */

    gp_data_read (1, 0, gp, m, 0);
  }

  /* Allocate space for hyperparameters, latent values, latent covariance matrix,
     and noise variances. */

  h->total_hypers = gp_hyper_count(gp,m);
  h->hyper_block = chk_alloc (h->total_hypers, sizeof (double));

  gp_hyper_pointers (h, gp, m);

  if (gen_latent)
  { latent_values = chk_alloc (gp->N_outputs*N_train, sizeof(double));
    latent_cov = chk_alloc (N_train*N_train, sizeof(double));
    scratch = chk_alloc (N_train, sizeof(double));
  }

  if (gen_noise_var)
  { noise_variances = chk_alloc (data_spec->N_targets*N_train, sizeof(double));
  }

  /* Read last records in log file to see where to start, and to get random
     number state left after last network was generated. */

  index = log_gobble_last(&logf,&logg);

  if (logg.last_index<0) 
  { index = 0;
  }

  if (index>max_index)
  { fprintf(stderr,"Records up to %d already exist in log file\n",max_index);
    exit(1);
  }

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  /* Generate new records with each index and write them to the log file. */

  for ( ; index<=max_index; index++)
  {
    /* Generate hyperparameter values. */

    gp_prior_generate (h, gp, m, fix, scale_value, relevance_value);

    /* Generate noise variances. */

    if (gen_noise_var)
    { double a, p, n;
      int i, j;
      a = m->noise.alpha[2];
      for (j = 0; j<data_spec->N_targets; j++)
      { n = exp (2 * *h->noise[j]);
        for (i = 0; i<N_train; i++)
        { p = rand_gamma(a/2) / (a*n/2);
          noise_variances [data_spec->N_targets*i + j] = 1/p;
        }
      }
    }

    /* Generate latent values. */

    if (gen_latent)
    { int i, j;
      for (j = 0; j<gp->N_outputs; j++)
      { gp_train_cov (gp, m, h, j, noise_variances, latent_cov, 0, 0);
        if (!gp->has_jitter && !(gp->lin.flags[j]&Flag_drop))
        { for (i = 0; i<N_train; i++)
          { latent_cov[i*N_train+i] += 1e-6;
          }
        }
        if (!cholesky(latent_cov,N_train,0))
        { fprintf (stderr,
            "Problem computing Cholesky decomposition of covariance matrix\n");
          exit(1);
        }
        for (i = 0; i<N_train; i++) 
        { scratch[i] = rand_gaussian();
        }
        for (i = 0; i<N_train; i++)
        { latent_values[i*gp->N_outputs+j] = 
            inner_product (scratch, 1, latent_cov+i*N_train, 1, i+1);
        }
      }
    }

    /* Write out records, including random number state record if we've
       been using random numbers. */

    logf.header.type = 'S';
    logf.header.index = index;
    logf.header.size = h->total_hypers * sizeof (double);
    log_file_append (&logf, h->hyper_block);

    if (gen_latent)
    { logf.header.type = 'F';
      logf.header.index = index;
      logf.header.size = gp->N_outputs*N_train * sizeof (double);
      log_file_append (&logf, latent_values);
    }
  
    if (gen_noise_var)
    { logf.header.type = 'N';
      logf.header.index = index;
      logf.header.size = data_spec->N_targets*N_train * sizeof (double);
      log_file_append (&logf, noise_variances);
    }
  
    if (!fix || gen_latent || gen_noise_var)
    { logf.header.type = 'r';
      logf.header.index = index;
      logf.header.size = sizeof (rand_state);
      log_file_append (&logf, rand_get_state());
    }
  }

  log_file_close(&logf);

  exit(0);
}
