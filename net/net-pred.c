/* NET-PRED.C - Make predictions for for test cases using neural networks. */

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
#include "data.h"
#include "net.h"
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
static net_flags *flgs;
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

  if (op_p && m==0)
  { fprintf(stderr,"Illegal combination of options with data model\n");
    exit(1);
  } 

  if (op_l)
  { fprintf(stderr,"Option l is not applicable to network models\n");
    exit(1);
  }

  a  = logg.data['A'];
  p  = logg.data['P'];

  flgs = logg.data['F'];

  net_check_specs_present(a,p,0,m,sv);

  s->total_sigmas = net_setup_sigma_count(a,flgs,m);
  w->total_params = net_setup_param_count(a,flgs);

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
  
  net_setup_sigma_pointers (s, a, flgs, m);
  net_setup_param_pointers (w, a, flgs);

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

    if (m==0 || m->type!='V' || sv->hazard_type=='C')
    {
      net_func (&test_values[i], 0, a, flgs, w);
    }
    
    if (have_targets) 
    { 
      if (op_p && m!=0 && m->type=='V' && sv->hazard_type!='C')
      {
        double ot, ft, t0, t1, lp;
        int censored;
        int x;
        
        if (test_targets[i]<0)
        { censored = 1;
          ot = -test_targets[i];
        }
        else
        { censored = 0;
          ot = test_targets[i];
        }

        t0 = 0;
        t1 = sv->time[0];
        test_values[i].i[0] = sv->log_time ? log(t1) : t1;

        x = 0;
        test_log_prob[i] = 0;

        for (;;)
        {
          net_func (&test_values[i], 0, a, flgs, w);
          
          ft = ot>t1 ? -(t1-t0) : censored ? -(ot-t0) : (ot-t0);

          net_model_prob(&test_values[i], &ft,
                         &lp, 0, a, m, sv, s, 0);

          test_log_prob[i] += lp;

          if (ot<=t1) break;
 
          t0 = t1;
          x += 1;
          
          if (sv->time[x]==0) 
          { t1 = ot;
            test_values[i].i[0] = sv->log_time ? log(t0) : t0;
          }
          else
          { t1 = sv->time[x];
            test_values[i].i[0] = sv->log_time ? (log(t0)+log(t1))/2
                                                  : (t0+t1)/2;
          }
        }
      }
      else
      {
        net_model_prob (&test_values[i], 
                        test_targets + data_spec->N_targets * i, 
                        &test_log_prob[i],
                        0, a, m, sv, s, 0);
      }

      if (op_r)
      { for (j = 0; j<data_spec->N_targets; j++)
        { tr = &data_spec->trans[data_spec->N_inputs+j];
          test_log_prob[i] += log(tr->scale);
          if (tr->take_log)
          { test_log_prob[i] -= log (data_inv_trans 
                               (test_targets[data_spec->N_targets*i+j], *tr));
          }
        }
      }
    }
    
    net_model_guess (&test_values[i], test_targ_pred + i*M_targets, 
                     a, flgs, m, sv, w, s, 0);

    if (op_D)
    { net_model_guess (&test_values[i], test_targ_med + i*M_targets, 
                       a, flgs, m, sv, w, s, 2);
    }

    for (j = 0; j<M_targets; j++)
    { 
      if (op_r) 
      { 
        tr = &data_spec->trans[data_spec->N_inputs+j];

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

    if (op_d || op_q || op_Q)
    {
      for (k = 0; k<Median_sample; k++)
      {
        net_model_guess (&test_values[i], curr_targets, 
                         a, flgs, m, sv, w, s, 1);

        for (j = 0; j<M_targets; j++)
        { if (op_r)
          { curr_targets[j] = data_inv_trans (curr_targets[j], 
                                 data_spec->trans[data_spec->N_inputs+j]);
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
