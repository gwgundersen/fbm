/* MIX-PRED.C - Make predictions for test cases using mixture models. */

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
#include "rand.h"
#include "log.h"
#include "prior.h"
#include "model.h"
#include "data.h"
#include "numin.h"
#include "mix.h"
#include "mix-data.h"
#include "mc.h"
#include "pred.h"


/* NAME OF THIS MODULE. */

char *pred_app_name = "mix";


/* ADJUSTABLE CONSTANTS. */

#define N_new_params 10		/* Number of random parameter values to use to
				   approximate predictions from new component */

/* LOCAL VARIABLES. */

static mix_spec *mx;		/* Diffusion tree model specification */

static int N_targets;		/* Number of target values in dataset */

static mix_hypers *hyp;		/* Hyperparameters for mixture model */

static int have_indicators;	/* Are component indicators present? */
static int have_offsets;	/* Are component offsets present? */
static int have_noise_SD;	/* Are component noise variables present? */

static short *indicators;	/* Pointer to array of indicator variables */
static double *offsets;		/* Pointer to array of offset parameters */
static double *noise_SD;	/* Pointer to array of noise std. dev. */

static int Max_active;		/* Maximum number of active components */
static int N_active;		/* Number of currently active components */

static int *freq;		/* Pointer to array of component frequencies */
static double *pr;		/* Array of relative component probabilities */


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* SET SIZES FOR APPLICATION RECORDS. */

void pred_app_record_sizes (void)
{
  mix_record_sizes(&logg);
}


/* INITIALIZE APPLICATION PROCEDURES. */

void pred_app_init (void)
{
  int dt; 

  if (op_N)
  { fprintf(stderr,
      "Options for selecting components are not allowed for mix-pred\n"
     );
    exit(1);
  }

  if (m==0 || (m->type!='R' && m->type!='B'))
  { fprintf(stderr,
      "Mixtures are implemented only for real and binary data models\n");
    exit(1);
  }

  mx = logg.data['P'];

  mix_check_specs_present(mx,0,m);

  N_targets = mx->N_targets;

  data_spec = logg.data['D'];

  if (data_spec==0) 
  { fprintf(stderr,"Can't make predictions with no training data\n");
    exit(1);
  }

  if (logg.actual_size['D'] !=
         data_spec_size(data_spec->N_inputs,data_spec->N_targets))
  { fprintf(stderr,"Data specification record is the wrong size!\n");
    exit(1);
  }
  
  mix_data_read (1, 1, mx, m);
}


/* LOOK AT MIXTURE MODEL STORED AT THE CURRENT INDEX. Returns 1 if there really
   something here, zero if not. */

int pred_app_use_index (void)
{    
  double c, s;
  int i, j, t;

  /* See what we have sitting in the log file for this iteration. */

  hyp = logg.data['S'];

  have_indicators = logg.data['I']!=0 && logg.index['I']==logg.last_index;
  have_offsets    = logg.data['O']!=0 && logg.index['O']==logg.last_index;
  have_noise_SD   = logg.data['N']!=0 && logg.index['N']==logg.last_index;

  /* See if there's really something here. */

  if (logg.index['S']!=logg.last_index)
  { 
    return 0;
  }

  /* Check that other stuff is here, and the right size. */

  if (!have_indicators || !have_offsets || m && m->type=='R' && !have_noise_SD)
  { fprintf(stderr,"Missing indicators, offsets, or noise records\n");
    exit(1);
  }

  if (logg.actual_size['I'] != N_train * sizeof *indicators)
  { fprintf(stderr,"Record of component indicators is wrong size!\n");
    exit(1);
  }

  Max_active = mx->N_components==0 ? N_train : mx->N_components;
  N_active = 0;

  indicators = logg.data['I'];

  for (i = 0; i<N_train; i++)
  { if (indicators[i]>=N_active) 
    { N_active = indicators[i] + 1;
    }
  }

  if (N_active>Max_active)
  { fprintf(stderr,"Garbled record of component indicators!\n");
    exit(1);
  }
 
  if (logg.actual_size['O'] != N_active*N_targets * sizeof *offsets)
  { fprintf(stderr,"Wrong size of record of component offsets!\n");
    exit(1);
  }

  offsets = logg.data['O'];
 
  if (have_noise_SD 
       && logg.actual_size['N'] != N_active*N_targets * sizeof *noise_SD)
  { fprintf(stderr, 
     "Wrong size of record of component noise standard deviations!\n");
    exit(1);
  }

  if (have_noise_SD) noise_SD = logg.data['N'];

  /* Allocate space for component frequencies, and initialize them. */

  freq = chk_alloc (N_active, sizeof *freq);
  mix_freq (indicators, N_train, freq, N_active);

  /* Compute probabilities for represented and new components (and allocate 
     space for them).*/

  pr = chk_alloc (N_active+1, sizeof *pr);

  if (mx->N_components==0)
  { c = 0;
  }
  else
  { c = hyp->con;
    if (mx->con_prior.scale) c /= mx->N_components;
  }

  for (i = 0; i<N_active; i++)
  { pr[i] = freq[i] + c;
  }

  pr[N_active] = mx->N_components==0 ? hyp->con : c*(mx->N_components-N_active);

  s = 0;
  for (i = 0; i<=N_active; i++) s += pr[i];
  for (i = 0; i<=N_active; i++) pr[i] /= s;

  /* Find probabilities for test cases.  Start by adding probbilites for 
     generating test cases from represented components. */

  for (j = 0; j<N_test; j++) 
  { 
    double lp, target;

    for (i = 0; i<N_active; i++)
    {
      lp = log(pr[i]);

      for (t = 0; t<N_targets; t++)
      {
        target = test_targets[N_targets*j+t];
  
        if (m->type=='R') /* Model for real data */
        { double d, v;
          d = target - offsets[N_targets*i+t];
          v = noise_SD[N_targets*i+t] * noise_SD[N_targets*i+t];
          lp += - (d*d) / (2*v) - log(2*M_PI*v)/2;
        }
        else if (m->type=='B') /* Model for binary data */
        { double l;
          l = offsets[N_targets*i+t];
          lp += target==0 ? -log(1+exp(l)) : -log(1+exp(-l));
        }
        else /* Type of model that isn't handled */
        { abort();
        }
      }

      test_log_prob[j] = i==0 ? lp : addlogs(test_log_prob[j],lp);
    }
  }

  /* Now add probabilities for generating test cases from a new component,
     approximating this with N_new_params parameter values picked at random. 

     Note that since the targets are independent given the parameters,
     and the parameters for different targets are independent given the
     hyperparameters, the average of products over targets can be replaced
     by the product of averages, which is more efficient in high dimensions. */

  if (pr[N_active]>0)
  { 
    double *noff, *nnsd;

    noff = chk_alloc (N_new_params * N_targets, sizeof *noff);
    if (m->type=='R') nnsd = chk_alloc (N_new_params * N_targets, sizeof *nnsd);

    for (i = 0; i<N_new_params; i++)
    { for (t = 0; t<N_targets; t++)
      { noff[N_targets*i+t] = hyp->mean[t] + hyp->SD[t] * rand_gaussian();
        if (m->type=='R') 
        { nnsd[N_targets*i+t] 
            = prior_pick_sigma (hyp->noise[t], m->noise.alpha[2]);
        }
      } 
    }

    for (j = 0; j<N_test; j++) 
    { 
      double lp, lpt, lp0, target;

      lp = log(pr[N_active]);

      for (t = 0; t<N_targets; t++)
      {
        target = test_targets[N_targets*j+t];
        
        for (i = 0; i<N_new_params; i++)
        {
          if (m->type=='R') /* Model for real data */
          { double d, v;
            d = target - noff[N_targets*i+t];
            v = nnsd[N_targets*i+t] * nnsd[N_targets*i+t];
            lp0 = - (d*d) / (2*v) - log(2*M_PI*v)/2;
          }
          else if (m->type=='B') /* Model for binary data */
          { double l;
            l = noff[N_targets*i+t];
            lp0 = target==0 ? -log(1+exp(l)) : -log(1+exp(-l));
          }
          else /* Type of model that isn't handled */
          { abort();
          }

          lpt = i==0 ? lp0 : addlogs(lpt,lp0);
        }
  
        lp += lpt - log((double)N_new_params);
      }

      test_log_prob[j] = addlogs(test_log_prob[j],lp);
    }

    if (m->type=='R') free(nnsd); 
    free(noff);
  }

  free(pr);
  free(freq);

  return 1;
}


/* CLEAN UP WHEN END OF LOG FILE IS REACHED. */

void pred_app_finish_file (void)
{
  logg.index['P'] = -1;  /* So they won't be mistaken for records from */
  logg.index['S'] = -1;  /* the new log file.                          */
  logg.index['I'] = -1;
  logg.index['O'] = -1;
  logg.index['N'] = -1;
}
