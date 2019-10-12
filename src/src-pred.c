/* SRC-PRED.C - Make predictions for test cases using source location model. */

/* Copyright (c) 1995-2007 by Radford M. Neal 
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
#include "phi.h"
#include "rand.h"
#include "log.h"
#include "mc.h"
#include "src.h"
#include "numin.h"
#include "data.h"
#include "src-data.h"
#include "prior.h"
#include "model.h"
#include "pred.h"


/* NAME OF THIS MODULE. */

char *pred_app_name = "src";


/* LOCAL VARIABLES. */

static src_spec *src;	/* Specification for source model */
static det_spec *det;	/* Specification for detector noise */
static flow_spec *flow;	/* Specification for flow model */

static src_params *params;

static double guessp, targ, prb, *prb_dist1, *prb_dist2;
static double *meanp, *varp; 


/* SET SIZES FOR APPLICATION RECORDS. */

void pred_app_record_sizes (void)
{
  src_record_sizes(&logg);
}


/* INITIALIZE APPLICATION PROCEDURES. */

void pred_app_init (void)
{
  static model_specification m;

  if (op_r)
  { fprintf (stderr, 
             "Prediction of raw target values is not currently supported\n");
    exit(1);
  }

  src = logg.data['S'];
  det = logg.data['T'];
  flow = logg.data['F'];

  src_check_specs_present(src,det,flow);
 
  src_data_read (0, 1);
  data_spec->N_inputs = 4;  /* Fudge this to make the "pred" skeleton work */

  meanp = chk_alloc (N_test, sizeof(double));
  varp  = chk_alloc (N_test, sizeof(double));

  /* Make up fake model specification to keep the "pred" skeleton happy. */

  logg.data['M'] = &m;
  m.type = 'R';
}


/* LOOK AT SOURCE MODEL STORED AT THE CURRENT INDEX. Returns 1 if there really
   is something here, zero if not. */

int pred_app_use_index (void)
{    
  int want_variance;
  int i, j, k, q;

  /* See if there's really something here. */

  if (logg.index['q']!=logg.last_index)
  { 
    return 0;
  }

  /* See what we have sitting in the log file for this iteration. */

  if (logg.actual_size['q'] != src_params_length(src->highN))
  { fprintf(stderr,"Bad parameter record (%d %d %d)\n",
            logg.actual_size['q'],src_params_length(src->highN),
            src->highN);
    exit(1);
  }
  params = logg.data['q'];

  want_variance = !op_z && (op_p || op_d || op_q);

  /* Find the predictive mean/variance for the output in each test case. */

  for (i = 0; i<N_test; i++)
  { 
    meanp[i] = src_total (src, flow, params, test_inputs+4*i);

    if (want_variance)
    { varp[i] = exp(2*params->log_width);
      /* needs to be modified later for t-distributed noise */
    }
  }

  /* Use the predictive means and variances of outputs to make predictions
     for the targets in each test case. */

  for (i = 0; i<N_test; i++) 
  { 
    test_targ_pred[i] = meanp[i];

    if (op_D)
    { test_targ_med[i] = meanp[i];
    }

    if (have_targets && op_p) 
    { test_log_prob[i] = 
        log_phi ((test_targets[i]-meanp[i])/sqrt(varp[i])) - 0.5*log(varp[i]);
    }

    if (op_d || op_q)
    { for (k = 0; k<Median_sample; k++)
      { double v;
        v = varp[i];
        if (params->inv_df>0)
        { v *= (0.5/params->inv_df) / rand_gamma(0.5/params->inv_df);
        }
        guessp = meanp[i] + sqrt(v)*rand_gaussian();
        median_sample[i][0][ms_count+k] = guessp;
      }
    }
  }

  return 1;
}


/* CLEAN UP WHEN END OF LOG FILE IS REACHED. */

void pred_app_finish_file (void)
{
  logg.index['S'] = -1;  /* So they won't be mistaken for records from */
  logg.index['T'] = -1;  /* the new log file.                          */
  logg.index['F'] = -1;
  logg.index['q'] = -1;
}
