/* MIX-CASES.C - Generate cases from a mixture model. */

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
#include "rand.h"
#include "mix.h"
#include "mix-data.h"


static void usage(void)
{ fprintf (stderr,
    "Usage: mix-cases [ -h ] [ -p ] [ -i ] log-file index\n");
  fprintf (stderr,
    "                 output-file n-cases [ random-seed ]\n");
  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  mix_spec *mx;
  mix_hypers *h;
  model_specification *m;

  int N_targets;
  int N_components;

  log_file logf;
  log_gobbled logg;

  int show_hypers;
  int show_params;
  int show_indicators;

  int index, n_cases;

  char *output_file;
  FILE *output;

  int random_seed;

  int Max_active;
  int N_active;

  short *indicators;
  double *offsets;
  double *noise_SD;

  int *freq;
  double *pr;

  int l, i, n, t;

  /* Look at arguments. */

  show_hypers= 0;
  show_params = 0;
  show_indicators = 0;

  while (argc>1 && *argv[1]=='-')
  {
    if (strcmp(argv[1],"-h")==0)
    { show_hypers = 1;
    }
    else if (strcmp(argv[1],"-p")==0)
    { show_params = 1;
    }
    else if (strcmp(argv[1],"-i")==0)
    { show_indicators = 1;
    }
    else 
    { usage();
    }

    argv += 1;
    argc -= 1;
  }

  if (argc!=5 && argc!=6) usage();

  logf.file_name = argv[1];

  if ((index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0) 
  { usage();
  }

  output_file = argv[3];

  if ((n_cases = atoi(argv[4]))<=0) usage();

  random_seed = index;
  if (argc==6 && (random_seed = atoi(argv[5]))<=0) usage();

  /* Open log file and read specification. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  mix_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  mx = logg.data['P'];
  m = logg.data['M'];

  mix_check_specs_present(mx,1,m);  

  N_targets = mx->N_targets;
  N_components = mx->N_components;

  /* Read the desired records from the log file. */

  while (!logf.at_end && logf.header.index!=index)
  { log_file_forward(&logf);
  }

  if (logf.at_end)
  { fprintf(stderr,"That index does not appear in the log file\n");
    exit(1);
  }

  log_gobble(&logf,&logg);

  if (logg.index['S']!=index)
  { fprintf(stderr,"No hyperparametes stored with that index\n");
    exit(1);
  }

  h = logg.data['S'];

  /* Read training data, if any, and determine maximum number of components. */
  
  data_spec = logg.data['D'];

  if (data_spec!=0) 
  { mix_data_read (1, 0, mx, m);
  }
  else
  { N_train = 0;
  }

  Max_active = N_train + n_cases;
  if (N_components!=0 && N_components<Max_active)
  { Max_active = N_components;
  }

  /* Allocate space for indicators, look for them in the log file if we
     have training data, and compute component frequencies. */

  indicators = chk_alloc (N_train+n_cases, sizeof *indicators);

  if (N_train>0)
  { 
    if (logg.data['I']==0)
    { fprintf(stderr,"No component indictors stored in log file\n");
      exit(1);
    }

    if (logg.actual_size['I'] != N_train * sizeof *indicators)
    { fprintf(stderr,"Record of component indicators is wrong size!\n");
      exit(1);
    }

    for (i = 0; i<N_train; i++) 
    { indicators[i] = ((short*)logg.data['I'])[i];
    }
  }

  N_active = 0;

  for (i = 0; i<N_train; i++)
  { if (indicators[i]>=N_active) 
    { N_active = indicators[i] + 1;
    }
  }

  if (N_active>N_train || N_active>Max_active)
  { fprintf(stderr,"Garbled record of component indicators!\n");
    exit(1);
  }

  freq = chk_alloc (Max_active, sizeof *freq);
  pr = chk_alloc (Max_active+1, sizeof *pr);
  
  mix_freq (indicators, N_train, freq, N_active);

  /* Allocate arrays to hold component parameters, and look for them in
     the log file if we have training data. */

  offsets = chk_alloc (Max_active*N_targets, sizeof *offsets);

  noise_SD = m->type!='R' ? 0
           : chk_alloc (Max_active*N_targets, sizeof *noise_SD);

  if (N_train>0)
  {
    if (logg.data['O']==0) 
    { fprintf(stderr,"Missing offset parameters for components\n");
      exit(1);
    }

    if (logg.actual_size['O'] != N_active*N_targets * sizeof *offsets)
    { fprintf(stderr,"Record of offset parameters is the wrong size!\n");
      exit(1);
    }

    for (i = 0; i<N_active*N_targets; i++) 
    { offsets[i] = ((double*)logg.data['O'])[i];
    }
  }

  if (N_train>0 && m->type=='R')
  { 
    if (logg.data['N']==0) 
    { fprintf(stderr,
       "Missing noise standard deviation parameters for components\n");
      exit(1);
    }

    if (logg.actual_size['N'] != N_active*N_targets * sizeof *noise_SD)
    { fprintf(stderr,
       "Record of noise standard deviation parameters is wrong size!\n");
      exit(1);
    }

    for (i = 0; i<N_active*N_targets; i++) 
    { noise_SD[i] = ((double*)logg.data['N'])[i];
    }
  }

  /* Open the output file. */

  output = fopen(output_file,"w");
  if (output==NULL)
  { fprintf(stderr,"Can't create output file (%s)\n",output_file);
    exit(1);
  }

  /* Generate the new cases, writing them to the output file. */

  rand_seed (random_seed);

  for (l = N_train; l<N_train+n_cases; l++)
  {
    /* Compute probabilities for the various components, including one
       that is not active (unless we're at the limit), then pick a particular
       component for the new case. */

    for (i = 0; i<N_active; i++) 
    { pr[i] = freq[i];
      if (N_components!=0) pr[i] += h->con/N_components;
    }

    if (N_components==0 || N_active<N_components)
    { pr[N_active] = N_components==0 ? h->con 
                   : h->con * (N_components-N_active) / N_components;
      n = rand_pickd(pr,N_active+1);
    }
    else
    { n = rand_pickd(pr,N_active);
    }

    /* Generate parameters from the prior (given hyperparameters) if the
       component picked is a new one. */

    if (n==N_active) 
    {
      for (t = 0; t<N_targets; t++)
      { 
        offsets[N_targets*n+t] = h->mean[t] + h->SD[t]*rand_gaussian();

        if (m->type=='R')
        { noise_SD[N_targets*n+t] = 
            prior_pick_sigma(h->noise[t],m->noise.alpha[2]);
        }
      }

      freq[n] = 0;

      N_active += 1;
    }

    /* Record indicator for new case, and update component frequencies. */

    indicators[l] = n;
    freq[n] += 1;

    /* Generate target values based on the chosen component's parameters. */

    for (t = 0; t<N_targets; t++)
    { if (m->type=='B')
      { fprintf (output, " %d",
          rand_uniform() < 1/(1+exp(-offsets[N_targets*n+t])));
      }
      else if (m->type=='R')
      { fprintf (output, " %+.6e",
          offsets[N_targets*n+t] + noise_SD[N_targets*n+t]*rand_gaussian());
      }
      else 
      { abort();
      }
    }

    fprintf(output,"\n");
  }

  /* Print information on standard output, if asked to. */

  if (show_params || show_indicators)
  { mix_sort (indicators, N_train+n_cases, freq, N_active,
              offsets, noise_SD, N_targets);
  }

  if (show_hypers)
  { mix_print_hypers (mx, m, h);
    printf("\n");
  }

  if (show_params)
  { mix_print_params (mx, m, N_active, freq, offsets, noise_SD);
    printf("\n");
  }

  if (show_indicators)
  { mix_print_indicators (N_train+n_cases, indicators);
    printf("\n");
  }

  exit(0);
}
