/* NET-QUANTITIES.C - Module defining quantities for neural networks. */

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
#include "net.h"
#include "net-data.h"


/* NETWORK VARIABLES. */

static net_arch *arch;
static net_flags *flgs;
static model_specification *model;
static net_priors *priors;
static model_survival *surv;

static net_sigmas sigmas;
static net_params params;

static net_params grad;
static net_values deriv;

static int M_targets;
static double *target_guess;

static char test_inputs_spec[100];
static char test_targets_spec[100];

static int have_train_data;
static int have_test_data;
static int have_test_targets;


/* LOOK AT ARGUMENTS. */

void net_arguments
( char ***argvp
)
{ 
  if (**argvp!=0)
  { 
    strcpy(test_inputs_spec,**argvp);
    *argvp += 1;

    if (**argvp!=0)
    { strcpy(test_targets_spec,**argvp);
      *argvp += 1;
    }
  }
}


/* INITIALIZE AFTER FIRST RECORDS READ. */

void net_initialize
( log_gobbled *logg
)
{ 
  /* Check that required specification records are present. */

  arch   = logg->data['A'];
  flgs   = logg->data['F'];
  model  = logg->data['M'];
  priors = logg->data['P'];
  surv   = logg->data['V'];

  net_check_specs_present(arch,priors,0,model,surv);

  /* Check that network is present, and set up pointers. */

  sigmas.total_sigmas = net_setup_sigma_count(arch,flgs,model);
  params.total_params = net_setup_param_count(arch,flgs);

  sigmas.sigma_block = logg->data['S'];
  params.param_block = logg->data['W'];

  if (sigmas.sigma_block==0 || params.param_block==0 
   || logg->index['S']!=logg->index['W'])
  { fprintf(stderr,"Network stored in log file is apparently incomplete\n");
    exit(1);
  }

  if (logg->actual_size['S'] != sigmas.total_sigmas*sizeof(net_sigma)
   || logg->actual_size['W'] != params.total_params*sizeof(net_param))
  { fprintf(stderr,"Bad size for network record\n");
    exit(1);
  }

  net_setup_sigma_pointers (&sigmas, arch, flgs, model);
  net_setup_param_pointers (&params, arch, flgs);

  grad.total_params = params.total_params;
  grad.param_block = chk_alloc (grad.total_params, sizeof (net_param));
  net_setup_param_pointers (&grad, arch, flgs);

  net_setup_value_pointers (&deriv, 
    chk_alloc (net_setup_value_count(arch), sizeof (net_value)), arch);

  /* Read training and test data, if present. */

  have_train_data = have_test_data = have_test_targets = 0;
  data_spec = logg->data['D'];

  if (data_spec!=0)
  {
    have_train_data = 1;

    if (test_inputs_spec[0]!=0)
    { strcpy (data_spec->test_inputs, test_inputs_spec);
      strcpy (data_spec->test_targets, test_targets_spec);
    }
    
    have_test_data = data_spec->test_inputs[0]!=0;
    have_test_targets = data_spec->test_targets[0]!=0;

    net_data_free ();   
    net_data_read (1, 1, arch, model, surv);

    M_targets = arch->N_outputs;

    target_guess = chk_alloc (M_targets, sizeof *target_guess);
  }
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE. */

void net_available 
( quantities_described qd,
  log_gobbled *logg
)
{ 
  char model_type = model ? model->type : 0;

  char letter;
  int mod;

  int v, o, n, s;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    mod = qd[v].modifier;

    if (letter && qd[v].available==0 && (strchr("iItT",letter)==0 || mod!=-1))
    {
      if (strchr("ixoygtzlbav",letter)!=0 && !have_train_data
       || strchr("IXOYGTZLBAV",letter)!=0 && !have_test_data
       || strchr("ZBA",letter)!=0       && !have_test_targets)
      { qd[v].available = -1;
        continue;
      }

# if 0 /* Disable these */
      if (strchr("uU",letter)!=0)
      { qd[v].available = 1;
      }
# endif
      if (strchr("xXiI",letter)!=0)
      { qd[v].available = mod<arch->N_inputs ? 1 : -1; 
      }
      else if (strchr("oOyYgGzZtT",letter)!=0)
      { qd[v].available = mod<arch->N_outputs ? 1 : -1; 
      }
      else if (strchr("bBaA",letter)!=0)
      { qd[v].available = mod<arch->N_outputs ? 1 : -1;
      }
      else if (letter=='c' || letter=='C')
      { qd[v].available = model_type=='C' && qd[v].low==-1 ? 1 : -1;
      }
      else if (strchr("PlL",letter)!=0)
      { qd[v].available = mod==-1 ? 1 : -1;
      }
      else if (letter=='h')
      { qd[v].available = 
           net_setup_hyper_group(arch,flgs,mod,&o,&n,&s) ? 1 : -1;
        if (qd[v].available==1 && qd[v].low==-1 && s)
        { qd[v].available = -1;
        }
        if (qd[v].available==1 && qd[v].low!=-1)
        { if (qd[v].high==-1) qd[v].high = n-2;
          if (n<=1 || qd[v].high>n-2 || qd[v].low>n-2) qd[v].available = -1;
        }
      }
      else if (letter=='w')
      { qd[v].available = qd[v].low!=-1 && 
           net_setup_param_group(arch,flgs,mod,&o,&n,&s) ? 1 : -1;
        if (qd[v].available==1)
        { if (qd[v].high==-1) qd[v].high = n-1;
          if (qd[v].high>n-1 || qd[v].low>n-1) qd[v].available = -1;
        }
      }
      else if (letter=='n' || letter=='N')
      { qd[v].available = (model_type=='R' || model_type=='C')
                             && mod==-1 ? 1 : -1;
      }
      else if (letter=='v' || letter=='V')
      { qd[v].available = mod<arch->N_layers ? 1 : -1;
      }
      else if (letter=='W')
      { qd[v].available = 
           net_setup_param_group(arch,flgs,mod,&o,&n,&s) ? 1 : -1;
        if (qd[v].available==1 && qd[v].low!=-1)
        { if (qd[v].high==-1) qd[v].high = s-1;
          if (s==0 || qd[v].high>s-1 || qd[v].low>s-1) qd[v].available = -1;
        }
      }
      else if (letter=='M')
      { qd[v].available = arch->N_layers>0 && mod<arch->N_outputs ? 1 : -1;
      }

      if (qd[v].available<0) continue;

      if (strchr("ixoygtzlba",letter)!=0)
      { if (strchr("ixoytzg",letter)!=0 && qd[v].low==-1 || qd[v].high>=N_train)
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = N_train-1;
      }

      else if (strchr("IXOYGTZLBA",letter)!=0)
      { if (strchr("IXOYTZG",letter)!=0 && qd[v].low==-1 || qd[v].high>=N_test) 
        { qd[v].available = -1;
          continue;
        }
        if (strchr("TZLAB",letter)!=0 && !have_test_targets)
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = N_test-1;
      }

      else if (letter=='P')
      { if (qd[v].low!=-1)
        { qd[v].available = -1;
          continue;
        }
      }
 
      else if (letter=='n' || letter=='N')
      { if (model_type=='R')
        { if (qd[v].high>=arch->N_outputs)
          { qd[v].available = -1;
            continue;
          }
          if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = arch->N_outputs-1;
        }
        else
        { if (qd[v].low!=-1) 
          { qd[v].available = -1;
            continue;
          }
        }
      }
 
      else if (letter=='v' || letter=='V')
      { if (qd[v].high >= arch->N_hidden [mod<0 ? 0 : mod])
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1) 
        { qd[v].high = arch->N_hidden [mod<0 ?0 : mod] - 1;
        }
      }
   
      else if (letter=='M')
      { if (qd[v].high >= arch->N_hidden[arch->N_layers-1]) 
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1)
        { qd[v].high = arch->N_hidden[arch->N_layers-1];
        }
        if (qd[v].low==-1)
        { qd[v].low = qd[v].high = 0;
        }
      }

      else if (strchr("oO",letter)!=0)
      { if (model_type=='V' && surv->hazard_type!='C') 
        { qd[v].available = -1;
          continue;
        }
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE. */

void net_evaluate 
( quantities_described qd, 
  quantities_held *qh,
  log_gobbled *logg
)
{ 
  char model_type = model ? model->type : 0;

  int mod, low, high;
  net_values *cases;
  double *targets;
  int N_cases;
  int ev_train, ev_test;
  char letter;
  int v, i, j;

  if (logg->data['S']==0 || logg->index['S']!=logg->last_index
   || logg->data['W']==0 || logg->index['W']!=logg->last_index) 
  { fprintf(stderr,"  records missing\n"); return;
  }

  ev_train = ev_test = 0;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    low  = qd[v].low;
    high = qd[v].high;
    mod  = qd[v].modifier;

    if (letter && !qh->updated[v] && (strchr("iItT",letter)==0 || mod!=-1))
    {
      if (letter>='a' && letter<='z')
      { cases = train_values;
        targets = train_targets;
        N_cases = N_train;
      }
      else
      { cases = test_values;
        targets = test_targets;
        N_cases = N_test;
      }

      if (!ev_train 
      && (strchr("uoyglbacv",letter)!=0 || letter=='n' && model_type=='C'))
      { if (!have_train_data) abort();
        for (i = 0; i<N_train; i++)
        { net_func (&train_values[i], 0, arch, flgs, &params);
        }
        ev_train = 1;
      }

      if (!ev_test 
      && (strchr("OYGLBACV",letter)!=0 || letter=='N' && model_type=='C'))
      { if (!have_test_data) abort();
        for (i = 0; i<N_test; i++)
        { net_func (&test_values[i], 0, arch, flgs, &params);
        }
        ev_test = 1;
      }

      if (mod<0) mod = 0;

      switch (letter)
      { 
        case 'x': case 'X': case 'i': case 'I':
        { for (i = low; i<=high; i++)
          { qh->value[v][i-low] = cases[i].i[mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'o': case 'O':
        { for (i = low; i<=high; i++)
          { qh->value[v][i-low] = cases[i].o[mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'y': case 'Y':
        { for (i = low; i<=high; i++)
          { net_model_guess (&cases[i], target_guess, arch, flgs, model, surv,
                             &params, &sigmas, 0);
            qh->value[v][i-low] = target_guess[mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'g': case 'G':                                
        { for (i = low; i<=high; i++)
          { net_model_guess (&cases[i], target_guess, arch, flgs, model, surv,
                             &params, &sigmas, 1);
            qh->value[v][i-low] = target_guess[mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'z': case 'Z': case 't': case 'T':
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
        { net_prior_prob (&params, &sigmas, qh->value[v], 0, 
                          arch, flgs, priors, 0);
          qh->updated[v] = 1;
          break;
        }

        case 'l': case 'L':
        { 
          double p, p1;

          if (low==-1) p = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            if (model_type=='V' && surv->hazard_type=='P')
            { 
              double ot, ft, t0, t1, lp;
              int censored;
              int w;

              if (targets[i]<0)
              { censored = 1;
                ot = -targets[i];
              }
              else
              { censored = 0;
                ot = targets[i];
              }

              p1 = 0;

              t0 = 0;
              t1 = surv->time[0];
              cases[i].i[0] = surv->log_time ? log(t1) : t1;

              w = 0;

              for (;;)
              {
                net_func (&cases[i], 0, arch, flgs, &params);
          
                ft = ot>t1 ? -(t1-t0) : censored ? -(ot-t0) : (ot-t0);

                net_model_prob(&cases[i], &ft, &lp, 0, arch, model, surv, 
                               &sigmas, 2);
                p1 += lp;

                if (ot<=t1) break;
 
                t0 = t1;
                w += 1;
          
                if (surv->time[w]==0) 
                { t1 = ot;
                  cases[i].i[0] = surv->log_time ? log(t0) : t0;
                }
                else
                { t1 = surv->time[w];
                  cases[i].i[0] = surv->log_time ? (log(t0)+log(t1))/2
                                                 : (t0+t1)/2;
                }
              }

            }
            else
            { net_model_prob (&cases[i], targets + data_spec->N_targets*i,
                              &p1, 0, arch, model, surv, &sigmas, 0);
            }

            if (low==-1)
            { p += p1;
            }
            else
            { qh->value[v][i-low] = -p1;
            }
          }
 
          if (low==-1)
          { *qh->value[v] = - p / N_cases;
          }

          qh->updated[v] = 1;

          break;
        }

        case 'a': case 'A':
        { 
          double d, e, tv;
          int m;

          if (low==-1) e = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            net_model_guess (&cases[i], target_guess, arch, flgs, model, surv,
                             &params, &sigmas, 0);

            if (low!=-1) 
            { qh->value[v][i-low] = 0;
            }

            for (m = mod; m<=(qd[v].modifier==-1 ? M_targets-1 : mod); m++)
            { tv = model_type=='C' ? targets[i]==m : targets[M_targets*i+m];
              if (model_type=='V' && tv<0) tv = -tv;
              d = target_guess[m] - tv;
              if (low==-1) 
              { e += (d>0 ? d : -d);
              }
              else 
              { qh->value[v][i-low] += (d>0 ? d : -d);
              }
            }
          }
 
          if (low==-1)
          { *qh->value[v] = e / N_cases;
          }

          qh->updated[v] = 1;

          break;
        }

        case 'b': case 'B':
        { 
          double d, e, tv;
          int m;

          if (low==-1) e = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            net_model_guess (&cases[i], target_guess, arch, flgs, model, surv,
                             &params, &sigmas, 0);

            if (low!=-1) 
            { qh->value[v][i-low] = 0;
            }

            for (m = mod; m<=(qd[v].modifier==-1 ? M_targets-1 : mod); m++)
            { tv = model_type=='C' ? targets[i]==m : targets[M_targets*i+m];
              if (model_type=='V' && tv<0) tv = -tv;
              d = target_guess[m] - tv;
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

        case 'c': case 'C':
        { 
          double e, p;
          int c;

          e = 0;

          for (i = 0; i<N_cases; i++)
          { 
            net_model_guess (&cases[i], target_guess, arch, flgs, model, surv,
                             &params, &sigmas, 0);

            p = qd[v].modifier==-1 ? target_guess[0] 
                                   : target_guess[0]+target_guess[1];
            for (c = 1 + (qd[v].modifier!=-1); c<arch->N_outputs; c++)
            { if (target_guess[c]>p) p = target_guess[c];
            }
 
            e += 1-p;
          }

          *qh->value[v] = e / N_cases; 
          qh->updated[v] = 1;

          break;
        }

# if 0  /* Disable the "U" and "u" quantities - probably not generally useful */

        case 'u':
        { 
          double u;
          int o, n, s;
          int i, j;

          if (mod<1) 
          { o = 0; 
            n = params.total_params;
          }
          else
          { if (!net_setup_param_group (arch, flgs, mod, &o, &n, &s)) abort();
          }
          
          for (j = o; j<o+n; j++) 
          { grad.param_block[j] = 0;
          }

          for (i = 0; i<N_train; i++)
          { 
            net_model_prob (&train_values[i], 
                            train_targets + data_spec->N_targets*i,
                            0, &deriv, arch, model, surv, &sigmas, 2);
  
            net_back (&train_values[i], &deriv, arch->has_ti ? -1 : 0,
                      arch, flgs, &params);

            net_grad (&grad, &params, &train_values[i], &deriv, arch,flgs);
          }

          u = 0;
          
          for (j = o; j<o+n; j++) 
          { u += grad.param_block[j] * grad.param_block[j];
          }

          *qh->value[v] = u;
          qh->updated[v] = 1;
 
          break;
        }

        case 'U':
        { 
          double u;
          int o, n, s;
          int i, j;

          if (mod<1) 
          { o = 0; 
            n = params.total_params;
          }
          else
          { if (!net_setup_param_group (arch, flgs, mod, &o, &n, &s)) abort();
          }
          
          u = 0;
          
          for (i = 0; i<N_train; i++)
          { 
            for (j = o; j<o+n; j++) 
            { grad.param_block[j] = 0;
            }

            net_model_prob (&train_values[i], 
                            train_targets + data_spec->N_targets*i,
                            0, &deriv, arch, model, surv, &sigmas, 2);
  
            net_back (&train_values[i], &deriv, arch->has_ti ? -1 : 0,
                      arch, flgs, &params);

            net_grad (&grad, &params, &train_values[i], &deriv, arch, flgs);

            for (j = o; j<o+n; j++) 
            { u += grad.param_block[j] * grad.param_block[j];
            }
          }

          *qh->value[v] = u;
          qh->updated[v] = 1;
 
          break;
        }

# endif

        case 'h':
        { 
          int o, n, s;

          if (!net_setup_hyper_group (arch, flgs, mod, &o, &n, &s)) abort();

          if (low==-1) 
          { *qh->value[v] = sigmas.sigma_block[o]; 
          }
          else
          { for (i = low; i<=high; i++)
            { qh->value[v][i-low] = sigmas.sigma_block[o+i+1];
            }
          }
 
          qh->updated[v] = 1;

          break;
        }

        case 'w':
        { 
          int o, n, s;

          if (!net_setup_param_group (arch, flgs, mod, &o, &n, &s)) abort();

          for (i = low; i<=high; i++)
          { qh->value[v][i-low] = params.param_block[o+i];
          }
 
          qh->updated[v] = 1;
        
          break;
        }

        case 'W':
        { 
          int o, n, s, j;
          double t, q;

          if (!net_setup_param_group (arch, flgs, mod, &o, &n, &s)) abort();

          if (low==-1)
          { q = 0;
            for (j = 0; j<n; j++)
            { t = params.param_block[o+j];
              q += t*t;
            }
            *qh->value[v] = sqrt(q/n);
          }
          else
          { for (i = low; i<=high; i++)
            { q = 0;
              for (j = 0; j<(n/s); j++)
              { t = params.param_block[o+i*(n/s)+j];
                q += t*t;
              }
              qh->value[v][i-low] = sqrt(q/(n/s));
            }
          }
 
          qh->updated[v] = 1;
        
          break;
        }

        case 'n': case 'N':
        { 
          if (model_type=='R')
          {
            if (low==-1)
            { *qh->value[v] = letter=='n' ? *sigmas.noise_cm 
                            : *sigmas.noise_cm * *sigmas.noise_cm;
            }
            else
            { for (i = low; i<=high; i++)
              { qh->value[v][i-low] = letter=='n' ? sigmas.noise[i]
                                    :  sigmas.noise[i] *  sigmas.noise[i];
              }
            }
            qh->updated[v] = 1;
            break;
          }
          else if (model_type=='C')
          {
            double e;
            int c;
  
            e = 0;
  
            for (i = 0; i<N_cases; i++)
            { 
              net_model_guess (&cases[i], target_guess, arch, flgs, model, surv,
                               &params, &sigmas, 0);
  
              for (c = 0; c<arch->N_outputs; c++)
              { e -= target_guess[c]==0 ? 0 
                      : target_guess[c] * log(target_guess[c]);
              }
            }
  
            *qh->value[v] = e / N_cases; 
            qh->updated[v] = 1;
  
            break;
          }
        }

        case 'v': case 'V':
        { 
          double tv, m, q, h;
 
          if (low==-1) tv = 0;

          for (j = (low==-1 ? 0 : low); 
               j <= (low==-1 ? arch->N_hidden[mod]-1 : high); 
               j++)
          { 
            m = q = 0;

            for (i = 0; i<N_cases; i++)
            { h = cases[i].h[mod][j];
              m += h;
              q += h*h;
            }

            q /= N_cases;
            m /= N_cases;

            if (low==-1) 
            { tv += q-m*m;
            }
            else 
            { qh->value[v][j-low] = sqrt(q-m*m);
            }
          }

          if (low==-1) 
          { *qh->value[v] = sqrt(tv) / arch->N_hidden[mod];
          }

          qh->updated[v] = 1;
          break;
        }

        case 'M':
        {
          double sorted[1000], t;
          int i, j;

          for (i = 0; i<arch->N_hidden[arch->N_layers-1]; i++)
          { t = params.ho[arch->N_layers-1][i*arch->N_outputs+mod];
            if (t<0) t = -t;
            for (j = i; j>0 && t>sorted[j-1]; j--)
            { sorted[j] = sorted[j-1];
            }
            sorted[j] = t;
          }

          for (i = low; i<=high; i++)
          { qh->value[v][i-low] = sorted[i];
          }
 
          qh->updated[v] = 1;
 
          break;
        }
      }
    }
  }
}
