/* MIX-GEN.C - Program to generate mixture model hyperparameters. */

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
#include "rand.h"


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  static mix_hypers hypers;  /* Static so initialized to zero, just in case */

  mix_spec *mx;

  model_specification *m;

  int max_index, index, fixarg, fix;
  double SD_value;

  log_file logf;
  log_gobbled logg;

  int t;
 
  /* Look at program arguments. */

  max_index = 0;
  SD_value = 0;

  fixarg = argc>2 && strcmp(argv[2],"fix")==0 ? 2
         : argc>3 && strcmp(argv[3],"fix")==0 ? 3
         : argc;

  fix = fixarg!=argc;

  if (argc<2 || argc!=fixarg && argc!=fixarg+1 && argc!=fixarg+2
   || fixarg>2 && (max_index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0
   || fixarg+1<argc && (SD_value = atof(argv[fixarg+1]))<=0)
  { fprintf(stderr,
    "Usage: mix-gen log-file [ max-index ] [ \"fix\" [ SD-value ] ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 1);

  log_gobble_init(&logg,0);
  mix_record_sizes(&logg);
  logg.req_size['r'] = sizeof (rand_state);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  mx = logg.data['P'];
  m = logg.data['M'];

  mix_check_specs_present(mx,0,m);

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

    if (fix)
    { 
      mix_hyper_init (mx, m, &hypers);

      if (SD_value!=0)
      { 
        if (mx->SD_prior.alpha[0]!=0)
        { hypers.SD_cm = SD_value;
        }

        for (t = 0; t<mx->N_targets; t++) 
        { hypers.SD[t] = mx->SD_prior.alpha[1]!=0 ? SD_value : hypers.SD_cm;
        }

        for (t = 0; t<mx->N_targets; t++)
        { hypers.mean[t] = 0;
        }

        if (m!=0 && m->type=='R')
        {
          if (m->noise.alpha[0]!=0)
          { hypers.noise_cm = SD_value;
          }

          for (t = 0; t<mx->N_targets; t++)
          { hypers.noise[t] = m->noise.alpha[1]!=0 ? SD_value : hypers.noise_cm;
          }
        }
      }

    }
    else /* Not fixed */
    {
      hypers.con = prior_pick_sigma(mx->con_prior.width,mx->con_prior.alpha[0]);

      hypers.SD_cm = 
         prior_pick_sigma (prior_width_scaled(&mx->SD_prior,mx->N_targets),
                           mx->SD_prior.alpha[0]);

      for (t = 0; t<mx->N_targets; t++) 
      { hypers.SD[t] = prior_pick_sigma(hypers.SD_cm,mx->SD_prior.alpha[1]);
      }

      for (t = 0; t<mx->N_targets; t++)
      { hypers.mean[t] = mx->mean_prior.width * rand_gaussian();
      }

      if (m!=0 && m->type=='R')
      {
        hypers.noise_cm = prior_pick_sigma(m->noise.width,m->noise.alpha[0]);

        for (t = 0; t<mx->N_targets; t++)
        { hypers.noise[t] = prior_pick_sigma(hypers.noise_cm,m->noise.alpha[1]);
        }
      }
      
    }

    /* Write out records, including random number state record if we've
       been using random numbers. */

    logf.header.type = 'S';
    logf.header.index = index;
    logf.header.size = sizeof hypers;
    log_file_append (&logf, &hypers);

    if (!fix)
    { logf.header.type = 'r';
      logf.header.index = index;
      logf.header.size = sizeof (rand_state);
      log_file_append (&logf, rand_get_state());
    }
  }

  log_file_close(&logf);

  exit(0);
}
