/* NET-REJ.C - Generate networks from posterior by rejection sampling. */

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
#include "net.h"
#include "net-data.h"
#include "rand.h"


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  net_arch *a;
  net_flags *flgs;
  net_priors *p;
  model_specification *m;
  model_survival *sv;

  net_sigmas  sigmas, *s = &sigmas;
  net_params  params, *w = &params;

  int max_index, limit, rejects, tries, index;

  log_file logf;
  log_gobbled logg;

  double log_prob, lp;

  int i;

  /* Look at program arguments. */

  limit = -1;

  if (argc!=3 && argc!=4
   || (max_index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0
   || argc>3 && (limit = atoi(argv[3]))<=0)
  { fprintf(stderr,"Usage: net-rej log-file max-index [ rejection-limit ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read network architecture and priors. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  net_record_sizes(&logg);
  logg.req_size['r'] = sizeof (rand_state);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }
  
  a = logg.data['A'];
  p = logg.data['P'];
  m = logg.data['M'];
  sv = logg.data['V'];

  flgs = logg.data['F'];

  net_check_specs_present(a,p,1,m,sv);

  /* Allocate space for parameters and hyperparameters. */

  s->total_sigmas = net_setup_sigma_count(a,flgs,m);
  w->total_params = net_setup_param_count(a,flgs);

  s->sigma_block = chk_alloc (s->total_sigmas, sizeof (net_sigma));
  w->param_block = chk_alloc (w->total_params, sizeof (net_param));

  net_setup_sigma_pointers (s, a, flgs, m);
  net_setup_param_pointers (w, a, flgs);

  /* Read last records in log file to see where to start, and to get random
     number state left after last network was generated. */

  index = log_gobble_last(&logf,&logg);

  if (logg.last_index<0) 
  { index = 0;
  }

  if (index>max_index)
  { fprintf(stderr,"Networks up to %d already exist in log file\n",max_index);
    exit(1);
  }

  if (logg.data['r']!=0) 
  { rand_use_state(logg.data['r']);
  }

  /* Read training data. */

  if ((data_spec = logg.data['D'])==0)
  { fprintf(stderr,"No data specification given\n");
    exit(1);
  }

  net_data_read (1, 0, a, m, sv);

  /* Generate networks from the posterior and write them to the log file. */

  tries = 0;
  rejects = 0;

  for ( ; index<=max_index; index++)
  {
    for (;;)
    {
      net_prior_generate (w, s, a, flgs, m, p, 0, 0, 0);

      if (m->type=='R') 
      { 
        double alpha;
 
        alpha = m->noise.alpha[0];

        *s->noise_cm = alpha==0 ? m->noise.width 
          : prior_pick_sigma (m->noise.width 
                             * sqrt(alpha/(N_train*a->N_outputs+alpha)),
                            N_train*a->N_outputs + alpha);

        alpha = m->noise.alpha[1];

        for (i = 0; i<a->N_outputs; i++)
        { s->noise[i] = alpha==0 ? *s->noise_cm 
           : prior_pick_sigma (*s->noise_cm * sqrt(alpha/(N_train+alpha)), 
                             N_train + alpha);
        }
      }

      log_prob = 0; 

      for (i = 0; i<N_train; i++)
      {
        if (m->type=='V' && sv->hazard_type!='C')
        { fprintf(stderr,
"Rejection sampling for non-constant survival models is not implemented yet\n");
          exit(1);
        }
        else
        { 
          net_func (&train_values[i], 0, a, flgs, w);
          net_model_prob (&train_values[i], 
                          train_targets + data_spec->N_targets*i,
                          &lp, 0, a, m, sv, s, 2);
          log_prob += lp;
        }
      }

      tries += 1;

      if (log_prob>0)
      { fprintf(stderr,
          "ERROR: Acceptance probability is greater than one (%.3le)\n",
          exp(log_prob));
        exit(1);
      }

      if (rand_uniform()<exp(log_prob)) break;

      rejects += 1;
    
      if (rejects==limit)
      { fprintf(stderr,"Rejection limit reached\n");
        goto done;
      }
    } 

    logf.header.type = 'S';
    logf.header.index = index;
    logf.header.size = s->total_sigmas * sizeof (net_sigma);
    log_file_append (&logf, s->sigma_block);

    logf.header.type = 'W';
    logf.header.index = index;
    logf.header.size = w->total_params * sizeof (net_param);
    log_file_append (&logf, w->param_block);

    logf.header.type = 'r';
    logf.header.index = index;
    logf.header.size = sizeof (rand_state);
    log_file_append (&logf, rand_get_state());
  }

  log_file_close(&logf);

done:

  fprintf (stderr, "Acceptance rate: %.1le\n", (double)(tries-rejects)/tries);

  exit(0);
}
