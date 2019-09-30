/* DFT-QUANTITIES.C - Module defining quantities for diffusion trees. */

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
#include "dft.h"
#include "dft-data.h"


/* DIFFUSION TREE MODEL SPECIFICATIONS. */

static dft_spec *dft;
static model_specification *model;

static int have_train_data;


/* FLAGS TO HELP ALGORITHMS. */

static int *flags;	/* Flags, indexed by numbers from 0 to N_train */
static int flags_n;	/* Number of training cases flags is allocated for */


/* PROCEDURES. */

static int split_number (int, int, int *, int *);


/* INITIALIZE AFTER FIRST RECORDS READ. */

void dft_initialize
( log_gobbled *logg
)
{ 
  /* Check that required specification records are present. */

  dft    = logg->data['P'];
  model  = logg->data['M'];

  dft_check_specs_present(dft,0,model);

  /* Read training and test data, if present. */

  have_train_data = 0;
  data_spec = logg->data['D'];
  N_train = 0;

  if (data_spec!=0)
  {
    have_train_data = 1;

    dft_data_free ();   
    dft_data_read (1, 1, dft, model);
  }
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE. */

void dft_available 
( quantities_described qd,
  log_gobbled *logg
)
{ 
  char model_type = model ? model->type : 0;

  char letter;
  int low, high, mod;
  int dt, a[3];
  int v;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    low  = qd[v].low;
    high = qd[v].high;
    mod = qd[v].modifier;

    if (letter && qd[v].available==0 && (letter!='t' || mod!=-1))
    {
      switch (letter)
      {
        case 't':
        { if (mod==-1) break;
          qd[v].available = -1;
          if (mod>=1 && mod<=dft->N_targets 
           && low<=N_train && high<=N_train)
          { if (high==-1) qd[v].high = N_train;
            if (low<=0) qd[v].low = 1;
            qd[v].available = 1;
          }
          break;
        }

        case 'o':
        { qd[v].available = -1;
          if ((mod==-1 || mod>=0 && mod<=dft->N_targets) 
           && low<=N_train && high<=N_train)
          { if (high==-1) qd[v].high = N_train;
            if (low<=0) qd[v].low = 1;
            if (mod==-1) qd[v].modifier = 1;
            qd[v].available = 1;
          }
          break;
        }

        case 'v': case 'V':
        { qd[v].available = -1;
          if ((mod==-1 || mod>=1 && mod<=dft->N_trees) 
            && low<=dft->N_targets && high<=dft->N_targets)
          { if (low!=-1 && high==-1) qd[v].high = dft->N_targets; 
            if (low==0) qd[v].low = 1;
            if (mod==-1) qd[v].modifier = 1;
            qd[v].available = 1;
          }
          break;
        }

        case 'c':
        { qd[v].available = -1;
          if ((mod==-1 || mod>=1 && mod<=dft->N_trees) && low>=0 && high<=2)
          { if (mod==-1) qd[v].modifier = 1;
            if (high==-1) qd[v].high = 2;
            qd[v].available = 1;
          }
          break;
        }

        case 'n': case 'N':
        { qd[v].available = -1;
          if (model_type=='R' && mod==-1 
           && low<=dft->N_targets && high<=dft->N_targets)
          { if (low!=-1 && high==-1) qd[v].high = dft->N_targets; 
            if (low==0) qd[v].low = 1;
            qd[v].available = 1;
          }
          break;
        }

        case 'd':
        { if (mod==-1) break;
          qd[v].available = -1;
          if (mod>=1 && mod<=dft->N_trees && high<=N_train-1)
          { if (low!=-1 && high==-1) qd[v].high = N_train-1;
            if (low!=-1) qd[v].low = 1;
            qd[v].available = 1;
          }
          break;
        }

        case 'D':
        { if (mod==-1) break;
          qd[v].available = -1;
          if (mod>=1 && mod<=dft->N_trees && low<=N_train && high<=N_train)
          { if (low!=-1 && high==-1) qd[v].high = N_train;
            if (low==0) qd[v].low = 1;
            qd[v].available = 1;
          }
          break;
        }

        case 'a':
        { qd[v].available = -1;
          if (low==-1 && split_number(mod,2,&dt,a)
           && a[0]>0 && a[0]<=N_train 
           && a[1]>0 && a[1]<=N_train)
          { qd[v].available = 1;
          }
          break;
        }

        case 'A':
        { qd[v].available = -1;
          if (low==-1 && split_number(mod,3,&dt,a)
           && a[0]>0 && a[0]<=N_train
           && a[1]>0 && a[1]<=N_train 
           && a[2]>0 && a[2]<=N_train
           && a[0]!=a[1] && a[0]!=a[2] && a[1]!=a[2])
          { qd[v].available = 1;
          }
          break;
        }

        default:
        { break;
        }
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE. */

void dft_evaluate 
( quantities_described qd, 
  quantities_held *qh,
  log_gobbled *logg
)
{ 
  char model_type = model ? model->type : 0;

  int mod, low, high;
  double *inputs;
  double *targets;
  dft_hypers *hyp;
  char letter;
  int v, i, j;
  int dt, a[3];

  int have_trees;		/* Are the trees present? */
  int have_latent;		/* Are latent vectors present? */
  int have_locations;		/* Are the node locations present? */

  int *parents;			/* Parents of nodes */
  double *divt;			/* Divergence times */
  double *latent;		/* Latent vectors */
  double *locations;		/* Node locations */

  dft_state st;			/* Pointers into state */

  if (logg->data['S']==0 || logg->index['S']!=logg->last_index)
  { fprintf(stderr,"  records missing\n"); return;
  }

  hyp = logg->data['S'];

  have_trees     = logg->data['R']!=0 && logg->index['R']==logg->last_index
                && logg->data['T']!=0 && logg->index['T']==logg->last_index;
  have_latent    = logg->data['L']!=0 && logg->index['L']==logg->last_index;
  have_locations = logg->data['N']!=0 && logg->index['N']==logg->last_index;

  parents   = logg->data['R'];
  divt      = logg->data['T'];
  latent    = logg->data['L'];
  locations = logg->data['N'];

  dft_setup_state (dft, model, st, hyp, parents, divt, locations, 0, N_train);

  if (flags_n!=N_train+1)
  { if (flags!=0) free(flags);
    flags = chk_alloc (N_train+1, sizeof *flags);
    flags_n = N_train+1;
  }

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    low  = qd[v].low;
    high = qd[v].high;
    mod  = qd[v].modifier;

    if (letter && !qh->updated[v] && (letter!='t' || mod!=-1))
    {
      inputs = train_inputs;
      targets = train_targets;

      switch (letter)
      { 
        case 't':
        { if (mod==-1) break;
          for (i = low; i<=high; i++)
          { qh->value[v][i-low] 
              = train_targets[data_spec->N_targets*(i-1)+mod-1];
          }
          qh->updated[v] = 1;
          break;
        }

        case 'o':
        { if (have_latent)
          { j = dft->N_targets*(low-1) + mod-1;
            for (i = low; i<=high; i++)
            { qh->value[v][i-low] = latent[j];
              j += dft->N_targets;
            }
            qh->updated[v] = 1;
          } 
          break;
        }

        case 'v': case 'V':
        { if (low==-1)
          { *qh->value[v] = hyp->d_SD_cm[mod-1];
            if (letter=='V')
            { *qh->value[v] *= *qh->value[v];
            }
          }
          else
          { for (i = low; i<=high; i++)
            { qh->value[v][i-low] = st[mod-1].d_SD[i-1];
              if (letter=='V')
              { qh->value[v][i-low] *= qh->value[v][i-low];
              }
            }
          }
          qh->updated[v] = 1;
          break;
        }

        case 'c':
        { for (i = low; i<=high; i++)
          { qh->value[v][i-low] = i==0 ? hyp->c0[mod-1] 
                                : i==1 ? hyp->c1[mod-1]
                                : hyp->c2[mod-1];
            qh->updated[v] = 1;
          }
          break;
        }

        case 'n': case 'N':
        { if (low==-1) 
          { *qh->value[v] = hyp->noise_cm;
            if (letter=='N')
            { *qh->value[v] *= *qh->value[v];
            }
          }
          else
          { for (i = low; i<=high; i++)
            { qh->value[v][i-low] = st->noise[i-1];
              if (letter=='N')
              { qh->value[v][i-low] *= qh->value[v][i-low];
              }
            }
          }
          qh->updated[v] = 1;
          break;
        }

        case 'd':
        { if (mod==-1) break;
          if (have_trees)
          { double *d;
            int n, m;
            d = divt + (N_train-1)*(mod-1);
            if (low==-1)
            { *qh->value[v] = 0;
              for (j = 0; j<N_train-1; j++)
              { *qh->value[v] += d[j];
              }
              *qh->value[v] /= N_train-1;
            }
            else
            { n = 0;
              for (j = 0; j<N_train-1; j++)
              { if (n<high-low+1 || d[j]<qh->value[v][n-1])
                { if (n<high-low+1) 
                  { n += 1;
                  }
                  for (m = n-1; m>0 && d[j]<qh->value[v][m-1]; m--)
                  { qh->value[v][m] = qh->value[v][m-1];
                  }
                  qh->value[v][m] = d[j];
                }
              }
              if (n<high-low+1) abort();
            }
            qh->updated[v] = 1;
          }
          break;
        }

        case 'D':
        { if (mod==-1) break;
          if (have_trees)
          { if (low==-1) 
            { *qh->value[v] = 0;
              for (i = 1; i<=N_train; i++)
              { int d;
                d = 0;
                j = i;
                while (st[mod-1].parents[j]!=0)
                { j = st[mod-1].parents[j];
                  if (j >= 0 || j < -(N_train-1)) abort();
                  d += 1;
                }
                *qh->value[v] += d;
              }
              *qh->value[v] /= N_train;
            }
            else
            { for (i = low; i<=high; i++)
              { int d;
                d = 0;
                j = i;
                while (st[mod-1].parents[j]!=0)
                { j = st[mod-1].parents[j];
                  if (j >= 0 || j < -(N_train-1)) abort();
                  d += 1;
                }
                qh->value[v][i-low] = d;
              }
            }
            qh->updated[v] = 1;
          }
          break;
        }

        case 'a':
        { if (!split_number(mod,2,&dt,a)) abort();
          if (have_trees)
          { int k;
            if (a[0]==a[1])
            { *qh->value[v] = 1;
            }
            else
            { dt -= 1;
              k = st[dt].parents[a[0]];
              while (k!=0)
              { flags[-k] = 0;
                k = st[dt].parents[k];
              }
              k = st[dt].parents[a[1]];
              while (k!=0)
              { flags[-k] = 1;
                k = st[dt].parents[k];
              }
              k = st[dt].parents[a[0]];
              while (flags[-k]==0)
              { k = st[dt].parents[k];
               if (k==0) abort();
              }
              *qh->value[v] = st[dt].divt[-k];
            }
            qh->updated[v] = 1;
          }
          break;
        }

        case 'A':
        { if (!split_number(mod,3,&dt,a)) abort();
          dt -= 1;
          if (have_trees)
          { int k;
            k = st[dt].parents[a[0]];
            while (k!=0)
            { flags[-k] = 0;
              k = st[dt].parents[k];
            }
            k = st[dt].parents[a[1]];
            while (k!=0)
            { flags[-k] = 1;
              k = st[dt].parents[k];
            }
            k = st[dt].parents[a[2]];
            while (k!=0)
            { flags[-k] = 2;
              k = st[dt].parents[k];
            }
            k = st[dt].parents[a[0]];
            while (flags[-k]==0)
            { k = st[dt].parents[k];
              if (k==0) abort();
            }
            *qh->value[v] = flags[-k]==1;
            qh->updated[v] = 1;
          }
          break;
        }
 
        default:
        { break;
        }
      }
    }
  }
}


/* SEPARATE NUMBER INTO PARTS.  Takes a number whose decimal representaton 
   is of one of the forms ta, tab, tabc, etc., where t is a one-digit number 
   (not zero), and a, b, c, etc. are numbers with the same number of digits,
   and returns the numerical values of t, a, b, c, etc.  The value of the
   function is 1 if the number is valid, 0 if it is not (eg, abc does not
   have a length that is a multiple of three). */

static int split_number
( int n,		/* Number to split */
  int k,		/* Number of parts in addition to the leading digit */
  int *t,		/* Place to store leading digit */
  int *p		/* Place to store remaining parts */
)
{ 
  static int pow10[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
  int s, i;

  if (n<=0)
  { return 0;
  }

  s = (int) (1e-30 + log((double)n)/log(10.0));

  if (s==0 || s%k!=0 || s/k>7) 
  { return 0;
  }

  for (i = 0; i<k; i++)
  { p[k-i-1] = n % pow10[s/k];
    n /= pow10[s/k];
  }

  *t = n;

  return 1;
}
