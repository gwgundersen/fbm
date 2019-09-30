/* NET-QUANTITIES.C - Module defining quantities for neural networks. */

/* Copyright (c) 1995 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "misc.h"
#include "log.h"
#include "quantities.h"
#include "net.h"
#include "data.h"
#include "net-data.h"


/* NETWORK VARIABLES. */

static net_arch *arch;
static net_priors *priors;

static net_sigmas sigmas;
static net_params params;

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
  int i, j;

  /* Check that required specification records are present. */

  if ((arch = logg->data['A'])==0)
  { fprintf(stderr,"No architecture specification in log file\n");
    exit(1);
  }

  if ((priors = logg->data['P'])==0)
  { fprintf(stderr,"No prior specification in log file\n");
    exit(1);
  }

  /* Check that network is present, and set up pointers. */

  sigmas.total_sigmas = net_setup_sigma_count(arch);
  params.total_params = net_setup_param_count(arch);

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

  net_setup_sigma_pointers (&sigmas, arch);
  net_setup_param_pointers (&params, arch);

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
    net_data_read (1, 1, arch);

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
  char letter;
  int mod, low;
  int v, o, n, s;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    mod = qd[v].modifier;

    if (letter && qd[v].available==0)
    {
      if (strchr("xoygzlbav",letter)!=0 && !have_train_data
       || strchr("XOYGZLBAV",letter)!=0 && !have_test_data
       || strchr("ZBA",letter)!=0      && !have_test_targets)
      { qd[v].available = -1;
        continue;
      }

      if (strchr("xX",letter)!=0)
      { qd[v].available = mod<arch->N_inputs ? 1 : -1; 
      }
      else if (strchr("oOyYgGzZ",letter)!=0)
      { qd[v].available = mod<arch->N_outputs ? 1 : -1; 
      }
      else if (strchr("bBaA",letter)!=0)
      { qd[v].available = arch->data_model!='C' && mod<arch->N_outputs ? 1 : -1;
      }
      else if (letter=='c' || letter=='C')
      { qd[v].available = arch->data_model=='C' && qd[v].low==-1 ? 1 : -1;
      }
      else if (strchr("PlL",letter)!=0)
      { qd[v].available = mod==-1 ? 1 : -1;
      }
      else if (letter=='h')
      { qd[v].available = net_setup_hyper_group(arch,mod,&o,&n,&s) ? 1 : -1;
        if (qd[v].available==1 && qd[v].low==-1 && s)
        { qd[v].available = -1;
        }
        if (qd[v].available==1 && qd[v].low!=-1)
        { if (qd[v].high==-1) qd[v].high = n-2;
          if (n<=1 || qd[v].high>n-2 || qd[v].low>n-2) qd[v].available = -1;
        }
      }
      else if (letter=='w')
      { qd[v].available = 
          qd[v].low!=-1 && net_setup_param_group(arch,mod,&o,&n,&s) ? 1 : -1;
        if (qd[v].available==1)
        { if (qd[v].high==-1) qd[v].high = n-1;
          if (qd[v].high>n-1 || qd[v].low>n-1) qd[v].available = -1;
        }
      }
      else if (letter=='n' || letter=='N')
      { qd[v].available = (arch->data_model=='R' || arch->data_model=='C')
                             && mod==-1 ? 1 : -1;
      }
      else if (letter=='v' || letter=='V')
      { qd[v].available = mod<arch->N_layers ? 1 : -1;
      }
      else if (letter=='W')
      { qd[v].available = net_setup_param_group(arch,mod,&o,&n,&s) ? 1 : -1;
        if (qd[v].available==1 && qd[v].low!=-1)
        { if (qd[v].high==-1) qd[v].high = s-1;
          if (s==0 || qd[v].high>s-1 || qd[v].low>s-1) qd[v].available = -1;
        }
      }
      else if (letter=='M')
      { qd[v].available = arch->N_layers>0 && mod<arch->N_outputs ? 1 : -1;
      }

      if (qd[v].available<0) continue;

      if (strchr("xoygzlba",letter)!=0)
      { if (strchr("xoyz",letter)!=0 && qd[v].low==-1 || qd[v].high>=N_train) 
        { qd[v].available = -1;
          continue;
        }
        if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = N_train-1;
      }

      else if (strchr("XOYGZLBA",letter)!=0)
      { if (strchr("XOYZ",letter)!=0 && qd[v].low==-1 || qd[v].high>=N_test) 
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
      { if (arch->data_model=='R')
        { if (qd[v].high>=arch->N_outputs)
          { qd[v].available = -1;
            continue;
          }
          if (qd[v].low!=-1 && qd[v].high==-1) qd[v].high = arch->N_outputs;
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

    if (letter && !qh->updated[v])
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
      && (strchr("oyglbacv",letter)!=0 || letter=='n' && arch->data_model=='C'))
      { if (!have_train_data) abort();
        for (i = 0; i<N_train; i++)
        { net_func (&train_values[i], 0, arch, &params);
        }
        ev_train = 1;
      }

      if (!ev_test 
      && (strchr("OYGLBACV",letter)!=0 || letter=='N' && arch->data_model=='C'))
      { if (!have_test_data) abort();
        for (i = 0; i<N_test; i++)
        { net_func (&test_values[i], 0, arch, &params);
        }
        ev_test = 1;
      }

      low  = qd[v].low;
      high = qd[v].high;
      mod  = qd[v].modifier;

      if (mod<0) mod = 0;

      switch (letter)
      { 
        case 'x': case 'X':
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
          { net_model_guess (&cases[i], target_guess, arch, priors, &sigmas, 0);
            qh->value[v][i-low] = target_guess[mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'g': case 'G':                                
        { for (i = low; i<=high; i++)
          { net_model_guess (&cases[i], target_guess, arch, priors, &sigmas, 1);
            qh->value[v][i-low] = target_guess[mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'z': case 'Z':
        { for (i = low; i<=high; i++)
          { qh->value[v][i-low] = 
              arch->data_model=='C' && qd[v].modifier>=0 ? targets[i]==mod
               : targets[data_spec->N_targets*i+mod];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'P':
        { net_prior_prob (&params, &sigmas, qh->value[v], 0, arch, priors, 0);
          qh->updated[v] = 1;
          break;
        }

        case 'l': case 'L':
        { 
          double p, p1;

          if (low==-1) p = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            net_model_prob (&cases[i], targets + data_spec->N_targets*i,
                            &p1, 0, arch, priors, &sigmas, 0);
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
          double d, e;
          int m;

          if (low==-1) e = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            net_model_guess (&cases[i], target_guess, arch, priors, &sigmas, 0);

            if (low!=-1) 
            { qh->value[v][i-low] = 0;
            }

            for (m = mod; m<=(qd[v].modifier==-1 ? M_targets-1 : mod); m++)
            { d = target_guess[m] - targets[M_targets*i+m];
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
          double d, e;
          int m;

          if (low==-1) e = 0;

          for (i = (low==-1 ? 0 : low); i <= (low==-1 ? N_cases-1 : high); i++)
          { 
            net_model_guess (&cases[i], target_guess, arch, priors, &sigmas, 0);

            if (low!=-1) 
            { qh->value[v][i-low] = 0;
            }

            for (m = mod; m<=(qd[v].modifier==-1 ? M_targets-1 : mod); m++)
            { d = target_guess[m] - targets[M_targets*i+m];
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
            net_model_guess (&cases[i], target_guess, arch, priors, &sigmas, 0);

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

        case 'h':
        { 
          int o, n, s;

          if (!net_setup_hyper_group (arch, mod, &o, &n, &s)) abort();

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

          if (!net_setup_param_group (arch, mod, &o, &n, &s)) abort();

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

          if (!net_setup_param_group (arch, mod, &o, &n, &s)) abort();

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
          if (arch->data_model=='R')
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
          else
          {
            double e;
            int c;
  
            e = 0;
  
            for (i = 0; i<N_cases; i++)
            { 
              net_model_guess (&cases[i], target_guess, arch, 
                               priors, &sigmas, 0);
  
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
