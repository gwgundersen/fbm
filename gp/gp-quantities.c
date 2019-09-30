/* GP-QUANTITIES.C - Module defining quantities for Gaussian processes. */

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
#include "quantities.h"
#include "prior.h"
#include "model.h"
#include "data.h"
#include "gp.h"
#include "gp-data.h"


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* GAUSSIAN PROCESS VARIABLES. */

static gp_spec *gp;
static model_specification *model;
static model_survival *surv;

static gp_hypers hypers;

static int M_targets;

static int have_train_data;


/* INITIALIZE AFTER FIRST RECORDS READ. */

void gp_initialize
( log_gobbled *logg
)
{ 
  /* Check that required specification records are present. */

  gp     = logg->data['P'];
  model  = logg->data['M'];
  surv   = logg->data['V'];

  gp_check_specs_present(gp,0,model,surv);

  /* Check that hyperparameters present; set up pointers. */

  hypers.total_hypers = gp_hyper_count(gp,model);
  hypers.hyper_block = logg->data['S'];

  if (hypers.hyper_block==0)
  { fprintf(stderr,"Gaussian process records in log file are incomplete\n");
    exit(1);
  }

  if (logg->actual_size['S'] != hypers.total_hypers*sizeof(double))
  { fprintf(stderr,"Bad size for Gaussian process record\n");
    exit(1);
  }

  gp_hyper_pointers (&hypers, gp, model);

  /* Read training and test data, if present. */

  have_train_data = 0;
  data_spec = logg->data['D'];

  if (data_spec!=0)
  {
    have_train_data = 1;

    gp_data_free ();   
    gp_data_read (1, 1, gp, model, surv);

    M_targets = gp->N_outputs;
  }
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE. */

void gp_available 
( quantities_described qd,
  log_gobbled *logg
)
{ 
  char model_type = model ? model->type : 0;

  char letter;
  int mod;
  int v;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    mod = qd[v].modifier;

    if (letter && qd[v].available==0 && (strchr("it",letter)==0 || mod!=-1))
    {
      if (strchr("ixoytzlbavV",letter)!=0 && !have_train_data)
      { qd[v].available = -1;
        continue;
      }

      if (letter=='x' || letter=='i')
      { qd[v].available = mod<gp->N_inputs ? 1 : -1; 
      }
      else if (strchr("oytz",letter)!=0)
      { qd[v].available = mod<gp->N_outputs ? 1 : -1; 
      }
      else if (strchr("ba",letter)!=0)
      { qd[v].available = mod<gp->N_outputs ? 1 : -1;
      }
      else if (letter=='P')
      { qd[v].available = mod==-1 ? 1 : -1;
      }
      else if (letter=='l')
      { qd[v].available = model_type!=0 && mod==-1 ? 1 : -1;
      }
      else if (letter=='S')
      { qd[v].available = 
          mod<=0 && gp->has_constant || mod>0 && mod<=gp->N_exp_parts ? 1 : -1;
      }
      else if (letter=='R')
      { qd[v].available = 
          mod<=0 && gp->has_linear || mod>0 && mod<=gp->N_exp_parts ? 1 : -1;
      }
      else if (letter=='M')
      { qd[v].available = qd[v].low!=-1 && mod<=gp->N_exp_parts ? 1 : -1;
      }
      else if (letter=='G')
      { qd[v].available = mod==-1 && gp->has_jitter ? 1 : -1;
      }
      else if (letter=='n' || letter=='N')
      { qd[v].available = model_type=='R' && mod==-1 ? 1 : -1;
      }
      else if (letter=='v' || letter=='V')
      { qd[v].available = model_type=='R' && mod<gp->N_outputs ? 1 : -1;
      }

      if (qd[v].available<0) continue;

      if (strchr("ixoytzlLbavV",letter)!=0)
      { if (strchr("ixoytzvV",letter)!=0 && qd[v].low==-1 
         || qd[v].high>=N_train) 
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = N_train-1;
      }

      else if (letter=='P' || letter=='S' || letter=='G')
      { if (qd[v].low!=-1)
        { qd[v].available = -1;
          continue;
        }
      }

      else if (letter=='R')
      { if (qd[v].low>=gp->N_inputs || qd[v].high>=gp->N_inputs)
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = gp->N_inputs-1;
      }

      else if (letter=='M')
      { if (qd[v].low>=gp->N_inputs || qd[v].high>=gp->N_inputs)
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = gp->N_inputs-1;
      }
 
      else if (letter=='n' || letter=='N')
      { if (qd[v].low>=gp->N_outputs || qd[v].high>=gp->N_outputs)
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = gp->N_outputs-1;
      }

      else if (letter=='o' || letter=='y')
      { if (model_type=='V' && surv->hazard_type!='C') 
        { qd[v].available = -1;
          continue;
        }
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE. */

void gp_evaluate 
( quantities_described qd, 
  quantities_held *qh,
  log_gobbled *logg
)
{ 
  char model_type = model ? model->type : 0;

  int mod, low, high;
  double *inputs;
  double *targets;
  double *latent_values;
  double *noise_variances;
  double pw;
  int N_cases;
  char letter;
  int v, i;

  if (logg->data['S']==0 || logg->index['S']!=logg->last_index)
  { fprintf(stderr,"  records missing\n"); return;
  }

  if (logg->data['S']!=hypers.hyper_block) abort();

  latent_values = logg->data['F']!=0 && logg->index['F']==logg->last_index
                   ? (double*) logg->data['F'] : 0;

  noise_variances = logg->data['N']!=0 && logg->index['N']==logg->last_index
                     ? (double*) logg->data['N'] : 0;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    low  = qd[v].low;
    high = qd[v].high;
    mod  = qd[v].modifier;

    if (letter && !qh->updated[v] && (strchr("it",letter)==0 || mod!=-1))
    {
      inputs = train_inputs;
      targets = train_targets;
      N_cases = N_train;

      if (mod<0) mod = 0;

      switch (letter)
      { 
        case 'x': case 'i':
        { for (i = low; i<=high; i++)
          { qh->value[v][i-low] = inputs[gp->N_inputs*i+mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'o':
        { if (latent_values!=0)
          { for (i = low; i<=high; i++)
            { qh->value[v][i-low] = latent_values[gp->N_outputs*i+mod];
            }
            qh->updated[v] = 1;
          }
          break;
        }

        case 'y':
        { 
          if (latent_values!=0)
          { 
            for (i = low; i<=high; i++)
            { 
              if (model_type=='B')
              { qh->value[v][i-low] = 
                  1 / (1+exp(-latent_values[gp->N_outputs*i+mod]));
              }
              else if (model_type=='N')
              { qh->value[v][i-low] = exp(latent_values[gp->N_outputs*i+mod]);
              }
              else if (model_type=='C')
              { double s;
                int m;
                s = 0;
                for (m = 0; m<gp->N_outputs; m++)
                { s += exp(latent_values[gp->N_outputs*i+m]);
                }
                qh->value[v][i-low] = exp(latent_values[gp->N_outputs*i+mod])/s;
              }
              else
              { qh->value[v][i-low] = latent_values[gp->N_outputs*i+mod];
              }
            }

            qh->updated[v] = 1;
          }
          break;
        }

        case 'z': case 't':
        { for (i = low; i<=high; i++)
          { qh->value[v][i-low] = 
              model_type=='C' && qd[v].modifier>=0 ? targets[i]==mod
               : targets[data_spec->N_targets*i+mod];
            if (model_type=='V' && qh->value[v][i-low]<0)
            { qh->value[v][i-low] = -qh->value[v][i-low];
            }
          }
          qh->updated[v] = 1;
          break;
        }

        case 'P':
        { *qh->value[v] = gp_log_prior(&hypers,gp,model,1);
          qh->updated[v] = 1;
          break;
        }

        case 'l':
        { 
          double l, e, s, tv, fv, nv;
          int m;

          if (latent_values==0 
           || noise_variances==0 && model_type=='R' && model->noise.alpha[2]!=0)
          { break;
          }

          if (low==-1) e = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            l = - gp_likelihood (&hypers, model, data_spec, 
              &targets[i*data_spec->N_targets], &latent_values[i*gp->N_outputs],
              noise_variances ? &noise_variances[i*data_spec->N_targets] : 0);
 
            if (low!=-1) 
            { qh->value[v][i-low] = l;
            }
            else
            { e += l;
            }
          }
 
          if (low==-1)
          { *qh->value[v] = e / N_cases;
          }

          qh->updated[v] = 1;

          break;
        }

        case 'a': 
        { 
          double d, e, s, tv, fv, gv;
          int m;

          if (latent_values==0) break;

          if (low==-1) e = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            if (low!=-1) 
            { qh->value[v][i-low] = 0;
            }

            if (model_type=='C')
            { s = 0;
              for (m = 0; m<gp->N_outputs; m++)
              { fv = latent_values[gp->N_outputs*i+m];
                s += exp(fv);
              }
            }

            for (m = mod; m<=(qd[v].modifier==-1 ? M_targets-1 : mod); m++)
            { tv = model_type=='C' ? targets[i]==m : targets[M_targets*i+m];
              if (model_type=='N' && tv<0) continue;
              fv = latent_values[M_targets*i+m];
              if (model_type=='R')
              { gv = fv;
              }
              else if (model_type=='B')
              { gv = 1/(1+exp(-fv));
              }
              else if (model_type=='N')
              { gv = exp(fv);
              }
              else if (model_type=='C')
              { gv = exp(fv)/s;
              }
              else
              { abort();
              }
              d = gv-tv;
              if (d<0) d = -d;
              if (low==-1)
              { e += d;
              }
              else 
              { qh->value[v][i-low] += d;
              }
            }
          }
 
          if (low==-1)
          { *qh->value[v] = e / N_cases;
          }

          qh->updated[v] = 1;

          break;
        }

        case 'b':
        { 
          double d, e, s, tv, fv, gv;
          int m;

          if (latent_values==0) break;

          if (low==-1) e = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            if (low!=-1) 
            { qh->value[v][i-low] = 0;
            }

            if (model_type=='C')
            { s = 0;
              for (m = 0; m<gp->N_outputs; m++)
              { fv = latent_values[gp->N_outputs*i+m];
                s += exp(fv);
              }
            }

            for (m = mod; m<=(qd[v].modifier==-1 ? M_targets-1 : mod); m++)
            { tv = model_type=='C' ? targets[i]==m : targets[M_targets*i+m];
              if (model_type=='N' && tv<0) continue;
              fv = latent_values[M_targets*i+m];
              if (model_type=='R')
              { gv = fv;
              }
              else if (model_type=='B')
              { gv = 1/(1+exp(-fv));
              }
              else if (model_type=='C')
              { gv = exp(fv)/s;
              }
              else if (model_type=='N')
              { gv = exp(fv);
              }
              else
              { abort();
              }
              d = gv-tv;
              if (low==-1)
              { e += d*d;
              }
              else 
              { qh->value[v][i-low] += d*d;
              }
            }
          }
 
          if (low==-1)
          { *qh->value[v] = e / N_cases;
          }

          qh->updated[v] = 1;

          break;
        }

        case 'S':
        {
          *qh->value[v] =
            exp(mod==0 ? *hypers.constant : *hypers.exp[mod-1].scale);

          qh->updated[v] = 1;
 
          break;
        }

        case 'R':
        {
          if (low==-1)
          { *qh->value[v] = 
              exp(mod==0 ? *hypers.linear_cm : *hypers.exp[mod-1].rel_cm);
          }
          else
          { for (i = low; i<=high; i++)
            { if (mod>0 && gp->exp[mod-1].flags[i]&Flag_omit
               || mod==0 && gp->linear_flags[i]&Flag_omit)
              { qh->value[v][i-low] = 0;
              }
              else
              { qh->value[v][i-low] = 
                  exp(mod==0 ? *hypers.linear[i] : *hypers.exp[mod-1].rel[i]);
              }
            }
          }

          qh->updated[v] = 1;
 
          break;
        }

        case 'M':
        {
          int m;

          for (i = low; i<=high; i++)
          {
            qh->value[v][i-low] = 
               qd[v].modifier>0 || !gp->has_linear 
                 || gp->linear_flags[i]&Flag_omit ? 0
               : 0.5 * exp(2 * *hypers.linear[i]);

            for (m = 1; m<=gp->N_exp_parts; m++)
            {
              if ((qd[v].modifier<0 || m==qd[v].modifier)
               && !(gp->exp[m-1].flags[i]&Flag_omit))
              { pw = gp->exp[m-1].power;
                qh->value[v][i-low] += exp(2 * *hypers.exp[m-1].scale)
                              * (1 - exp( - exp(pw * *hypers.exp[m-1].rel[i])));
              }
            }
          }

          qh->updated[v] = 1;
 
          break;
        }

        case 'G':
        {
          *qh->value[v] = exp(*hypers.jitter);

          qh->updated[v] = 1;

          break;
        }

        case 'n': case 'N':
        { 
          if (low==-1)
          { *qh->value[v] = exp(letter=='n' ? *hypers.noise_cm 
                                 : *hypers.noise_cm * 2);
          }
          else
          { for (i = low; i<=high; i++)
            { qh->value[v][i-low] = exp(letter=='n' ? *hypers.noise[i]
                                         : *hypers.noise[i] * 2);
            }
          }

          qh->updated[v] = 1;

          break;
        }
        case 'v': case 'V':
        { if (noise_variances!=0)
          { for (i = low; i<=high; i++)
            { qh->value[v][i-low] = 
                letter=='V' ? noise_variances[gp->N_outputs*i+mod]
                            : sqrt(noise_variances[gp->N_outputs*i+mod]);
            }
            qh->updated[v] = 1;
          }
          break;
        }
      }
    }
  }
}
