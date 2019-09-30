/* MIX-DISPLAY.C - Program to print parameters of a mixture model. */

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


static void usage(void)
{ fprintf (stderr,
    "Usage: mix-display [ -h ] [ -p ] [ -i ] log-file [ index ]\n");
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

  log_file logf;
  log_gobbled logg;

  int show_hypers;
  int show_params;
  int show_indicators;

  int index;

  short *indicators;
  double *offsets;
  double *noise_SD;
  int *freq;

  int N_active, N_cases;

  int i;

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

  if (show_hypers==0 && show_params==0 && show_indicators==0)
  { show_hypers = 1;
    show_params = 1;
  }

  index = -1;

  if (argc!=2 && argc!=3 
   || argc>2 && (index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0) 
  { usage();
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  mix_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  mx = logg.data['P'];
  m = logg.data['M'];
  
  if (mx==0)
  { fprintf(stderr,"No specification for mixture model in log file\n");
    exit(1);
  }

  /* Read the desired records from the log file. */

  if (index<0)
  { 
    log_gobble_last(&logf,&logg);

    if (logg.last_index<0)
    { fprintf(stderr,"No mixture model record in log file\n");
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

  if (logg.index['S']!=index)
  { fprintf(stderr,"No hyperparameters stored with that index\n");
    exit(1);
  }

  h = logg.data['S'];

  /* Print title. */

  printf("\nMIXTURE MODEL IN FILE \"%s\" WITH INDEX %d\n", 
         logf.file_name, index);

  /* Print things that were asked for. */

  if (show_hypers) 
  { mix_print_hypers (mx, m, h);
    printf("\n");
  }

  if (show_params)
  { offsets  = logg.index['O']==index ? logg.data['O'] : 0;
    noise_SD = logg.index['N']==index ? logg.data['N'] : 0;
    indicators = logg.index['I']==index ? logg.data['I'] : 0;
    if (indicators==0 || offsets==0 || m!=0 && m->type=='R' && noise_SD==0)
    { printf("\nNO RECORD OF MODEL PARAMETERS STORED AT THIS INDEX\n");
    }
    else
    { N_cases  = logg.actual_size['I'] / sizeof *indicators;
      N_active = N_cases==0 ? 0 : indicators[0]+1;
      for (i = 1; i<N_cases; i++) 
      { if (indicators[i]+1>N_active) 
        { N_active = indicators[i]+1;
        }
      }
      if (logg.actual_size['I'] != N_cases * sizeof *indicators
       || logg.actual_size['O'] != N_active*mx->N_targets * sizeof *offsets
       || m!=0 && m->type=='R' &&
           logg.actual_size['N'] != N_active*mx->N_targets * sizeof *noise_SD)
      { fprintf(stderr,"Records sizes are not right!\n");
        exit(1);
      }
      freq = chk_alloc (N_active, sizeof *freq);
      mix_freq (indicators, N_cases, freq, N_active);
      mix_print_params (mx, m, N_active, freq, offsets, noise_SD);
    }
    printf("\n");
  }

  if (show_indicators)
  { indicators = logg.index['I']==index ? logg.data['I'] : 0;
    if (indicators==0)
    { printf("\nNO COMPONENT INDICATORS STORED AT THIS INDEX\n");
    }
    else
    { N_cases = logg.actual_size['I'] / sizeof *indicators;
      mix_print_indicators (N_cases, indicators);
    }
    printf("\n");
  }

  exit(0);
}
