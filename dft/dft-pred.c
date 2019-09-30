/* DFT-PRED.C - Make predictions for test cases using diffusion trees. */

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
#include "phi.h"
#include "rand.h"
#include "log.h"
#include "prior.h"
#include "model.h"
#include "data.h"
#include "numin.h"
#include "dft.h"
#include "dft-data.h"
#include "mc.h"
#include "pred.h"


/* NAME OF THIS MODULE. */

char *pred_app_name = "dft";


/* LOCAL VARIABLES. */

static dft_spec *dft;		/* Diffusion tree model specification */

static int N_targets;		/* Number of target values in dataset */

static dft_hypers *hyp;		/* Hyperparameters for diffusion tree model */

static int have_trees;		/* Are the trees present? */
static int have_latent;		/* Are latent vectors present? */
static int have_locations;	/* Are the node locations present? */

static int *parents;		/* Parents of nodes */
static double *divt;		/* Divergence times */
static double *latent;		/* Latent vectors */
static double *locations;	/* Node locations */

static double *latent0;		/* Space allocated for above, for when they */
static double *locations0;	/*   don't come from the log file           */

static dft_tree_node *nodes;	/* Internal tree representation */

static double *terminals[Max_trees]; /* Terminal node locations (temporary) */

static dft_state st;		/* Pointers to components of state */

static double *test_mean;	/* Mean of latent vector for test case */
static double *test_var;	/* Variance of latent vector for test case */

static int *way;		/* Which way to go at branch points */


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* SET SIZES FOR APPLICATION RECORDS. */

void pred_app_record_sizes (void)
{
  dft_record_sizes(&logg);
}


/* INITIALIZE APPLICATION PROCEDURES. */

void pred_app_init (void)
{
  int dt; 

  if (op_N)
  { fprintf(stderr,
      "Options for selecting components are not allowed for dft-pred\n"
     );
    exit(1);
  }

  if (m!=0 && m->type!='R' && m->type!='B')
  { fprintf(stderr,
    "The data model used is not implemented for Dirichlet diffusion trees\n");
    exit(1);
  }

  dft = logg.data['P'];

  dft_check_specs_present(dft,0,m);

  N_targets = dft->N_targets;

  data_spec = logg.data['D'];

  if (data_spec==0) 
  { fprintf(stderr,"Can't make predictions with no training data\n");
    exit(1);
  }
  
  dft_data_read (1, 1, dft, m);

  for (dt = 0; dt<dft->N_trees; dt++)
  { terminals[dt] = chk_alloc (N_train*N_targets, sizeof (double));
  }

  test_mean = chk_alloc (N_targets, sizeof *test_mean);
  test_var  = chk_alloc (N_targets, sizeof *test_var);

  way = chk_alloc (N_train, sizeof *way);

  nodes = chk_alloc (dft->N_trees*N_train, sizeof *nodes);

  latent0    = chk_alloc(1,dft_latent_size(dft,N_train));
  locations0 = chk_alloc(1,dft_locations_size(dft,N_train));
}


/* LOOK AT DIFFUSION TREE STORED AT THE CURRENT INDEX. Returns 1 if there really
   something here, zero if not. */

int pred_app_use_index (void)
{    
  int dt, i, j, t;

  /* See what we have sitting in the log file for this iteration. */

  hyp = logg.data['S'];

  have_trees     = logg.data['T']!=0 && logg.index['T']==logg.last_index;
  have_latent    = logg.data['L']!=0 && logg.index['L']==logg.last_index;
  have_locations = logg.data['N']!=0 && logg.index['N']==logg.last_index;

  if (have_trees && (logg.data['R']==0 || logg.index['R']!=logg.last_index))
  { fprintf(stderr,"Missing record of parents!\n");
    exit(1);
  }

  if (!have_trees && (logg.data['R']==0 || logg.index['R']!=logg.last_index))
  { fprintf(stderr,"Missing record of divergence times!\n");
    exit(1);
  }

  if (hyp==0 && (have_trees || have_latent || have_locations))
  { fprintf(stderr,"Missing hyperparameter record!\n");
    exit(1);
  }

  if (!have_trees && (have_latent || have_locations))
  { fprintf(stderr,"Missing records describing tree structures!\n");
    exit(1);
  }

  dft_check_sizes(&logg,dft,m,N_train);

  /* See if there's really something here. */

  if (logg.index['S']!=logg.last_index || !have_trees)
  { 
    return 0;
  }

  /* Set up pointers to existing state, or to previously allocated space. */

  parents   = logg.data['R'];
  divt      = logg.data['T'];

  latent    = have_latent   ? logg.data['L'] : latent0;
  locations = have_locations? logg.data['N'] : locations0;

  /* Create nodes and set up state. */

  for (dt = 0; dt<dft->N_trees; dt++)
  { (void) dft_conv_tree (parents+dt*(2*N_train)+N_train-1, 
                          nodes+dt*N_train, N_train);
  }

  dft_setup_state (dft, m, st, hyp, parents, divt, 
                   locations, nodes, N_train);

  /* Generate latent vectors and locations for training cases, if not here. */

  if (!have_locations && N_train>0)
  { if (dft->N_trees>1)
    { fprintf(stderr,
      "Can't sample for absent locations when there's more than one tree\n");
      exit(1);
    }
    dft_gibbs_locations (dft, m, st, 1, have_latent ? latent : 0);
  }

  if (!have_latent && N_train>0)
  { dft_sample_latent (dft, m, st, 1, latent);
  }

  /* Generate locations for terminal nodes. */

  dft_sample_terminals (dft, st, latent, terminals);

  /* Simulate generation of N_train new cases, which follow paths towards
     each of the old cases (a stratification method).  Use these simulated
     cases to assign probabilities or densities to test cases. */

  for (i = 1; i<=N_train; i++)
  {
    /* Find mean and variance of latent vector for test cases, by summing
       simulations from all the trees. */

    for (t = 0; t<N_targets; t++)
    { test_mean[t] = 0; 
      test_var[t] = 0;
    }

    for (dt = 0; dt<dft->N_trees; dt++)
    { 
      double *vp, *vc, *vn, *v, *w;
      double s;
      int x, y;

      x = i;
      while (st[dt].parents[x]!=0)
      { y = st[dt].parents[x];
        if (chld(st[dt].nodes[-y],0)==x)
        { way[-y] = 0;
        }
        else if (chld(st[dt].nodes[-y],1)==x)
        { way[-y] = 1;
        }
        else 
        { abort();
        }
        x = y;
      }

      v = locations + (N_train-1)*N_targets*dt;
      w = terminals[dt];

      dft_gen_path (hyp, st, dt, x, &x, &s, way);
  
      vp = st[dt].parents[x]==0 ? 0 : v + (-st[dt].parents[x]-1)*N_targets;
      vc = x>0 ? w+(x-1)*N_targets : v+(-x-1)*N_targets;

      for (t = 0; t<N_targets; t++)
      { double am, bm, av, bv, sd;
        sd = st[dt].d_SD[t];
        am = vp==0 ? 0 : vp[t];
        av = (sd*sd) * (vp==0 ? s : s - st[dt].divt[-st[dt].parents[x]]);
        bm = vc[t];
        bv = (sd*sd) * (x>0 ? 1-s : st[dt].divt[-x] - s);
        if (av+bv<=0) /* just in case */
        { test_mean[t] += (am+bm)/2;
          test_var[t] += 0;
        }
        else
        { test_mean[t] += (am*bv+bm*av)/(av+bv);
          test_var[t] += (av*bv)/(av+bv);
        }
        test_var[t] += (1-s) * (sd*sd);

/*printf("i=%d dt=%d t=%d mean=%.4f var=%.4f\n",i,dt,t,test_mean[t],test_var[t]);*/

        if (test_var[t]<0) abort();
      }
    }

    /* Find probabilities for test cases based on means and variances for 
       latent values found above. */

    for (j = 0; j<N_test; j++) 
    { 
      double lp, target;

      lp = 0;

      for (t = 0; t<N_targets; t++)
      {
        target = test_targets[N_targets*j+t];

        if (m==0 || m->type=='R') /* Model for real data */
        { double d, v;
          d = target - test_mean[t];
          v = test_var[t];
          if (m!=0) v += st->noise[t]*st->noise[t];
          lp += - (d*d) / (2*v) - log(2*M_PI*v)/2;
/*printf("lp += %.2f\n",- (d*d) / (2*v) - log(2*M_PI*v)/2);*/
        }
        else if (m->type=='B') /* Model for binary data */
        { double l;
          l = test_mean[t] + sqrt(test_var[t]) * rand_gaussian();
          lp += target==0 ? -log(1+exp(l)) : -log(1+exp(-l));
        }
        else /* Type of model that isn't handled */
        { abort();
        }
      }

      test_log_prob[j] = i==1 ? lp : addlogs(test_log_prob[j],lp);
/*printf("partial %d %d: %.2f %.2f\n",i,j,lp,test_log_prob[j]);*/
    }
  }

  /* Divide by N_train to get final averaged probabilities. */

  for (j = 0; j<N_test; j++) 
  { test_log_prob[j] -= log((double)N_train);
/*printf("final %d: %.2f\n",j,test_log_prob[j]);*/
  }

  return 1;
}


/* CLEAN UP WHEN END OF LOG FILE IS REACHED. */

void pred_app_finish_file (void)
{
  logg.index['S'] = -1;  /* So they won't be mistaken for records from */
  logg.index['T'] = -1;  /* the new log file.                          */
  logg.index['R'] = -1;
  logg.index['L'] = -1;
  logg.index['N'] = -1;
}
