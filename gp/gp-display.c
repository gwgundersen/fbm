/* GP-DISPLAY.C - Program to print parameters of Gaussian process model. */

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
#include "gp.h"


static void usage(void)
{ fprintf(stderr,"Usage: gp-display [ -h ] [ -l ] [ -n ] log-file [ index ]\n");
  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  gp_spec *gp;
  model_specification *m;

  gp_hypers hypers, *h = &hypers;

  double *values, *variances;

  log_file logf;
  log_gobbled logg;

  int show_hypers;
  int show_values;
  int show_std_dev;

  int index;

  int l, i, n;

  /* Look at arguments. */

  show_hypers= 0;
  show_values = 0;
  show_std_dev = 0;

  while (argc>1 && *argv[1]=='-')
  {
    if (strcmp(argv[1],"-h")==0)
    { show_hypers = 1;
    }
    else if (strcmp(argv[1],"-l")==0)
    { show_values = 1;
    }
    else if (strcmp(argv[1],"-n")==0)
    { show_std_dev = 1;
    }
    else 
    { usage();
    }

    argv += 1;
    argc -= 1;
  }

  if (show_hypers==0 && show_values==0 && show_std_dev==0)
  { show_hypers = 1;
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
  gp_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  gp = logg.data['P'];
  m = logg.data['M'];
  
  if (gp==0)
  { fprintf(stderr,"No specification for Gaussian process in log file\n");
    exit(1);
  }

  h->total_hypers = gp_hyper_count(gp,m);

  logg.req_size['S'] = h->total_hypers * sizeof(double);

  /* Read the desired records from the log file. */

  if (index<0)
  { 
    log_gobble_last(&logf,&logg);

    if (logg.last_index<0)
    { fprintf(stderr,"No Gaussian process record in log file\n");
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

  h->hyper_block = logg.data['S'];

  gp_hyper_pointers (h, gp, m);

  /* Print title. */

  printf("\nGAUSSIAN PROCESS IN FILE \"%s\" WITH INDEX %d\n", 
         logf.file_name, index);

  /* Print values of the hyperparameters if asked. */

  if (show_hypers)
  {
    printf("\nHYPERPARAMETERS\n");

    if (gp->has_constant)
    { printf("\nConstant part:\n\n %10.3f\n",exp(*h->constant));
    }

    if (gp->has_linear)
    { printf("\nLinear part:\n\n %10.3f :",exp(*h->linear_cm));
      for (i = 0; i<gp->N_inputs; i++)   
      { if (i>0 && i%5==0)
        { printf("\n             ");
        }
        if (gp->linear_flags[i]&Flag_omit)
        { printf("       -   ");
        }
        else
        { printf(" %10.3f",exp(*h->linear[i]));
        }
      }
      printf("\n");
    }

    if (gp->has_jitter)
    { printf("\nJitter part:\n\n %10.3g\n",exp(*h->jitter));
    }

    if (gp->N_exp_parts>0)
    { printf("\nExponential parts:\n");
      for (l = 0; l<gp->N_exp_parts; l++)
      { printf ("\n %10.3f\n %10.3f :", exp(*h->exp[l].scale),
                                        exp(*h->exp[l].rel_cm));
        for (i = 0; i<gp->N_inputs; i++)
        { if (i>0 && i%5==0)
          { printf("\n             ");
          }
          if (gp->exp[l].flags[i]&Flag_omit)
          { printf("       -   ");
          }
          else
          { printf(" %10.3f",exp(*h->exp[l].rel[i]));
          }
        }
        printf("\n");
      }
    }

    if (m!=0 && m->type=='R')
    { printf("\nNoise levels:\n\n %10.3f :",exp(*h->noise_cm));
      for (i = 0; i<gp->N_outputs; i++)
      { if (i>0 && i%5==0)
        { printf("\n             ");
        }
        printf(" %10.3f",exp(*h->noise[i]));
      }
      printf("\n");
    }
  }


  /* Print latent values for training cases if asked. */

  if (show_values)
  {
    printf("\nLATENT VALUES FOR TRAINING CASES\n\n");
  
    if ((values = (double*) logg.data['F']) == 0
     || (n = logg.actual_size['F'] / (gp->N_outputs*sizeof(double))) == 0)
    { printf("No latent values stored\n");
    }
    else
    { for (i = 0; i<n; i++)
      { printf("%5d:",i);
        for (l = 0; l<gp->N_outputs; l++)
        { if (l>0 && l%10==0)
          { printf("\n      ");
          }
          printf(" %+6.2f",values[i*gp->N_outputs+l]);
        }
        printf("\n");
      }
    }
  }

  /* Print noise standard deviations for training cases if asked. */

  if (show_std_dev)
  {
    printf("\nNOISE STANDARD DEVIATIONS FOR TRAINING CASES\n\n");
  
    if ((variances = (double*) logg.data['N']) == 0
     || (n = logg.actual_size['N'] / (gp->N_outputs*sizeof(double))) == 0)
    { printf("No noise standard deviations stored\n");
    }
    else
    { for (i = 0; i<n; i++)
      { printf("%5d:",i);
        for (l = 0; l<gp->N_outputs; l++)
        { if (l>0 && l%10==0)
          { printf("\n      ");
          }
          printf(" %6.3f",sqrt(variances[i*gp->N_outputs+l]));
        }
        printf("\n");
      }
    }
  }

  printf("\n");

  exit(0);
}
