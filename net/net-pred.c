/* NET-PRED.C - Make predictions for for test cases using neural networks. */

/* Copyright (c) 1995, 1996, 1998, 1999 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 *
 * Modifications to allow selection of a given number of networks and to
 * allow specification of ranges by cpu-time are adapted from modifications
 * by Carl Edward Rasmussen, 1995.
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
#include "net.h"
#include "data.h"
#include "numin.h"
#include "net-data.h"
#include "mc.h"
#include "pred.h"


/* NAME OF THIS MODULE. */

char *pred_app_name = "net";


/* SHARED VARIABLE. */

double *test_inputs;


/* LOCAL VARIABLES. */

static net_arch *a;
static net_priors *p;
static net_sigmas sigmas, *s = &sigmas;
static net_params params, *w = &params;

static double *curr_targets;


/* SET SIZES FOR APPLICATION RECORDS. */

void pred_app_record_sizes (void)
{
  net_record_sizes(&logg);
}


/* INITIALIZE APPLICATION PROCEDURES. */

void pred_app_init (void)
{
  double *t;
  int i, j;

  if (op_l)
  { fprintf(stderr,"Option l is not applicable to network models\n");
    exit(1);
  }

  a  = logg.data['A'];
  p  = logg.data['P'];

  net_check_specs_present(a,p,0,m,sv);

  if (op_p && m!=0 && m->type=='V' && sv->hazard_type!='C')
  { fprintf(stderr,
"Option p is implemented for survival models only when the hazard is constant\n"
    );
    exit(1);
  }

  s->total_sigmas = net_setup_sigma_count(a,m);
  w->total_params = net_setup_param_count(a);

  logg.req_size['S'] = s->total_sigmas * sizeof(net_sigma);
  logg.req_size['W'] = w->total_params * sizeof(net_param);

  curr_targets = chk_alloc (M_targets, sizeof (double));

  net_data_read (0, 1, a, m, sv);

  test_inputs = chk_alloc (data_spec->N_inputs*N_test, sizeof (double));

  t = test_inputs;
  for (i = 0; i<N_test; i++)
  { for (j = m!=0 && m->type=='V' && sv->hazard_type!='C' ? 1 : 0; 
         j<a->N_inputs; j++)
    { *t++ = test_values[i].i[j];
    }
  }
}


/* LOOK AT NETWORK STORED AT THE CURRENT INDEX.  Returns 1 if there really
   is a network here, zero if not. */

int pred_app_use_index (void)
{    
  int i, j, k, l;

  if (logg.index['W']!=logg.last_index || logg.index['S']!=logg.last_index)
  { 
    return 0;
  }

  s->sigma_block = logg.data['S'];
  w->param_block = logg.data['W'];
  
  net_setup_sigma_pointers (s, a, m);
  net_setup_param_pointers (w, a);

  if (op_N)
  { for (l = 0; l<a->N_layers; l++)
    { if (a->has_ho[l] && !keep[l])
      { for (k = 0; k<a->N_hidden[l]*a->N_outputs; k++)
        { w->ho[l][k] = 0;
        }
      }
    }
    if (a->has_io)
    { for (k = 0; k<a->N_inputs*a->N_outputs; k++)
      { w->io[k] = 0;
      }
    }
    if (a->has_bo)
    { for (k = 0; k<a->N_outputs; k++)
      { w->bo[k] = 0;
      }
    }
  }
    
  for (i = 0; i<N_test; i++)
  { 
    net_func (&test_values[i], 0, a, w);
    
    if (have_targets) 
    { 
      net_model_prob (&test_values[i], 
                      test_targets + data_spec->N_targets * i, 
                      &test_log_prob[i],
                      0, a, m, sv, s, 0);
      if (op_r)
      { for (j = 0; j<data_spec->N_targets; j++)
        { tr = &data_spec->target_trans[j];
          test_log_prob[i] += log(tr->scale);
          if (tr->take_log)
          { test_log_prob[i] -= log (data_inv_trans 
                               (test_targets[data_spec->N_targets*i+j], *tr));
          }
        }
      }
    }
    
    net_model_guess (&test_values[i], test_targ_pred + i*M_targets, 
                     a, m, sv, w, s, 0);
    
    for (j = 0; j<M_targets; j++)
    { 
      if (op_r) 
      { 
        tr = &data_spec->target_trans[j];

        test_targ_pred[i*M_targets+j] = 
            data_inv_trans(test_targ_pred[i*M_targets+j], *tr);

        if (tr->take_log)
        {
          if (op_n && m!=0 && m->type=='R' && (m->noise.alpha[0]!=0 
               || m->noise.alpha[1]!=0 || m->noise.alpha[2]!=0))
          { fprintf(stderr,"Predictive mean is undefined\n");
            exit(1);
          }

          test_targ_pred[i*M_targets+j] *= exp (m->noise.width*m->noise.width
                                          / (tr->scale*tr->scale*2));
        }
      }
    }

    if (op_d || op_q)
    {
      rand_seed(101*logg.last_index+i);

      for (k = 0; k<Median_sample; k++)
      {
        net_model_guess (&test_values[i], curr_targets, 
                         a, m, sv, w, s, 1);

        for (j = 0; j<M_targets; j++)
        { if (op_r)
          { curr_targets[j] = data_inv_trans (curr_targets[j], 
                                        data_spec->target_trans[j]);
          }
          median_sample[i][j][ms_count+k] = curr_targets[j];
        }

      }
    }
  }

  return 1;
}


/* CLEAN UP WHEN END OF LOG FILE IS REACHED. */

void pred_app_finish_file (void)
{
  logg.index['W'] = -1;  /* So they won't be mistaken for records from */
  logg.index['S'] = -1;  /* the new log file.                          */
}
