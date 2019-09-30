/* DFT-PROB.C - Compute probabilities for diffusion tree models. */

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
#include "mc.h"
#include "data.h"
#include "prior.h"
#include "model.h"
#include "dft.h"
#include "dft-data.h"


/* TABLES OF PRECOMPUTED VALUES. */

static double *log_factorial;	/* Logs of factorial from 0 to N_train */
static double *sum_reciprocals;	/* Table of sum from i=1 to n of 1/i, for
				   n from 0 to N_train */

/* INITIALIZE THE TABLES OF PRECOMPUTED VALUES.*/

static void initialize_tables (void)
{
  int n;

  if (log_factorial==0) 
  { 
    log_factorial = chk_alloc (N_train+1, sizeof *log_factorial);

    log_factorial[0] = 0;

    for (n = 1; n<=N_train; n++)
    { log_factorial[n] = log_factorial[n-1] + log((double)n);
    }
  }

  if (sum_reciprocals==0)
  { 
    sum_reciprocals = chk_alloc (N_train+1, sizeof *sum_reciprocals);

    sum_reciprocals[0] = 0;
    for (n = 1; n<=N_train; n++)
    { sum_reciprocals[n] = sum_reciprocals[n-1] + 1.0/n;
    }
  }
}


/* FIND PROBABILITY OF DIVERGENCES.  Computes the log of the probability 
   that all paths in the given tree diverge from previous paths when they 
   do, and not before.  The probabilities for taking the correct branch are 
   not included.  The value can be -INFINITY if divergence doesn't happen
   by time 1 even though it should.

   For external use, the last two parameters should be zero.  (Internal 
   recursive calls set them to other values.) */

double dft_log_prob_div
( dft_spec *dft,		/* Diffusion tree model specification */
  dft_hypers *hyp,		/* Hyperparameters for diffusion tree model */
  int dt,			/* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  int root,			/* Root of sub-tree to look at, zero for all */
  double div0			/* Divergence time of branch to root */
)
{
  double lp, div1;
  int n;

  initialize_tables();

  if (root==0)
  { for (root = 1; st[dt].parents[root]!=0; root = st[dt].parents[root]) ;
  }

  if (root>0) return 0;

  n = totpts(st[dt].nodes[-root]);
  if (n<2) abort();

  div1 = st[dt].divt[-root];
  if (div1<div0) abort();

  if (div1>div0)
  { lp = - sum_reciprocals[n-1] * (dft_cdiv(hyp,dt,div1)-dft_cdiv(hyp,dt,div0));
  }

  lp += log (dft_div (hyp, dt, div1));

  lp += dft_log_prob_div (dft, hyp, dt, st, chld(st[dt].nodes[-root],0), div1);
  lp += dft_log_prob_div (dft, hyp, dt, st, chld(st[dt].nodes[-root],1), div1);

  return lp;
} 


/* FIND PROBABILITY OF PATHS TO SUBTREE.  Returns the log probability of
   the paths in the tree to the nodes in the subtree headed by x turning
   out the way they did, excluding the portions of these paths within 
   this subtree.  The probability is the product of the probabilities of 
   not diverging too soon and taking the right turns at each branch before 
   reaching x, and also of the probability density for the divergence of the 
   first path from previous paths.  The result may be -INFINITY, when
   the divergence time is exactly 1. */

double dft_log_prob_paths
( dft_spec *dft,		/* Diffusion tree model specification */
  dft_hypers *hyp,		/* Hyperparameters for diffusion tree model */
  int dt,			/* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  int x				/* Node heading subtree */
)
{ 
  double logprob, A;
  int y, nx, ny, nc;

  initialize_tables();

  if (st[dt].parents[x]>=0) abort();

  nx = x>0 ? 1 : totpts(st[dt].nodes[-x]);

  logprob = x>0 ? 0 : - sum_reciprocals[nx-1] *
                        (dft_cdiv(hyp,dt,st[dt].divt[-x]) 
                          - dft_cdiv(hyp,dt,st[dt].divt[-st[dt].parents[x]]));

  for (y = st[dt].parents[x]; y!=0; y = st[dt].parents[y])
  { 
    ny = totpts(st[dt].nodes[-y]);

    if (y==st[dt].parents[x])
    { logprob += log (dft_div(hyp,dt,st[dt].divt[-y]) / (ny-nx))
               + log_factorial[nx-1]
               - (log_factorial[ny-1] - log_factorial[ny-nx]);
    }
    else
    { logprob += (log_factorial[nc-1] - log_factorial[nc-nx-1])
               - (log_factorial[ny-1] - log_factorial[ny-nx-1]);
    }

    A = dft_cdiv (hyp, dt, st[dt].divt[-y]);
    if (st[dt].parents[y]!=0) 
    { A -= dft_cdiv (hyp, dt, st[dt].divt[-st[dt].parents[y]]);
    }

    logprob -= A * (sum_reciprocals[ny-1] - sum_reciprocals[ny-nx-1]);

    nc = ny;
  }

  if (0 && x>0) /* Check that can be enabled for debugging */
  { double logprob2;
    logprob2 = dft_log_prob_path(dft,hyp,dt,st,x);
    if (logprob2-logprob>0.00001 || logprob-logprob2>0.00001) 
    { abort();
    }
  }

  if (isnan(logprob)) abort();

  return logprob;
}


/* FIND PROBABILITY OF PATH TO TERMINAL NODE.  Returns the log probability 
   of the first path in the tree that goes through the node x.  The probability
   is the product of the probabilities of taking the right turns at each branch 
   before reaching x, and also of the probability density for the divergence 
   of the path to x from previous paths.  The result may be -INFINITY, when
   the divergence time is exactly 1.*/

double dft_log_prob_path
( dft_spec *dft,		/* Diffusion tree model specification */
  dft_hypers *hyp,		/* Hyperparameters for diffusion tree model */
  int dt,			/* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  int x				/* Node heading subtree */
)
{ 
  double logprob, A;
  int y, nx, ny, nc;

  if (st[dt].parents[x]>=0) abort();

  nx = x>0 ? 1 : totpts(st[dt].nodes[-x]);

  for (y = st[dt].parents[x]; y!=0; y = st[dt].parents[y])
  { 
    ny = totpts(st[dt].nodes[-y]) - nx;

    if (y==st[dt].parents[x])
    { logprob = log (dft_div(hyp,dt,st[dt].divt[-y]) / ny);
    }
    else
    { logprob += log ((double) nc / ny);
    }

    A = dft_cdiv (hyp, dt, st[dt].divt[-y]);
    if (st[dt].parents[y]!=0) 
    { A -= dft_cdiv (hyp, dt, st[dt].divt[-st[dt].parents[y]]);
    }

    logprob -= A / ny;

    nc = ny;
  }

  if (isnan(logprob)) abort();

  return logprob;
}


/* FIND PROBABILITY DENSITY FOR NODE POSITION.  Finds the probability density
   for the observed data along with the latent vectors and other node locations,
   if present, perhaps omitting factors that do not depend on the point where 
   the parent of x is located.  This probability density integrates over the 
   location of the parent of x, as well as over the locations of terminal 
   nodes, and of any other absent data.  

   If node locations are present, a new location for the parent of x is
   randomly generated. 

   If the final "lka" parameter is present, it holds likelihood functions for
   each node and each target variable, which are used rather than computing 
   likelihoods from scratch.  This is of significance only when there are no 
   node locations. 

   The log probability is returned in the "logprob" argument if the value is
   finite.  Infinite values result in "infprob" being set to a value greater
   than zero, giving the number of points where infinite density occurs.  The
   "logprob" argument will in this case be set to the value excluding these
   infinities. */

void dft_log_prob_node
( dft_spec *dft,		/* Diffusion tree model specification */
  int dt,			/* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  double *latent,		/* Latent vectors, or null */
  int x,			/* Node heading subtree */
  double inv_temp,		/* Inverse temperature */
  double *logprob,		/* Set to log probability */
  int *infprob,			/* Set to number of points of infinite density*/
  dft_likelihood *lka		/* Likelihoods for nodes, or null */
)
{
  int N_targets;
  int t;

  if (st[dt].parents[x]>=0) abort();

  N_targets = dft->N_targets;

  if (st[dt].locations!=0)
  {
    int px, gpx, sx, t;
    double t0, t1;

    px = st[dt].parents[x];
    gpx = st[dt].parents[px];
    sx = dft_sibling (st[dt].parents, st[dt].nodes, x);
    t0 = gpx==0 ? 0 : st[dt].divt[-gpx];
    t1 = st[dt].divt[-px];

    *logprob = 0;
    *infprob = 0;

    for (t = 0; t<N_targets; t++)
    { 
      double pmean, pvar, var, sd, postmean, postvar, lmean, lvar, gpv;

      /* Find prior mean and precision for the parent of x. */

      dft_node_likelihood (dft, dt, st, latent, sx, t, inv_temp, &lmean, &lvar);

      sd = st[dt].d_SD[t];
      
      pvar = t1-t0<=0 ? 0 : lvar * (sd*sd)*(t1-t0) / (lvar + (sd*sd)*(t1-t0));

      gpv = gpx==0 ? 0 : st[dt].locations[N_targets*(-gpx-1)+t];

      pmean = t1-t0<=0 ? gpv 
            : (lvar * gpv + (sd*sd)*(t1-t0) * lmean) / (lvar + (sd*sd)*(t1-t0));

      /* Find likelihood for the parent of x w.r.t. x. */

      dft_node_likelihood (dft, dt, st, latent, x, t, inv_temp, &lmean, &lvar);

      /* Compute log probability of observation. */

      var = lvar + pvar;

      *logprob -= log_sqrt_2pi;
      if (var>0) 
      { *logprob += - (pmean-lmean)*(pmean-lmean)/(2*var) - log(var)/2;
      }
      else
      { if (pmean==lmean)
        { *infprob += 1;
        }
        else
        { *logprob -= INFINITY;
        }
      }

      /* Pick a random location for the parent of x. */

      postvar = var==0 ? 0 : (pvar*lvar) / var;
      postmean = var==0 ? (pmean + lmean) / 2 : (pmean*lvar + lmean*pvar) / var;

      st[dt].locations[(-px-1)*N_targets + t] = 
                                       postmean + rand_gaussian()*sqrt(postvar);
    }
  }

  else /* node locations are not present */
  { 
    double mean, var, peak, v0;
    int root, infp;
  
    for (root = x; st[dt].parents[root]!=0; root = st[dt].parents[root]) ;

    *logprob = 0;
    *infprob = 0;

    for (t = 0; t<N_targets; t++)
    { if (lka==0)
      { dft_tree_likelihood (dft, dt, st, latent, t, root, inv_temp,
                             &mean, &var, &peak, &infp, 0);
      }
      else
      { dft_likelihood *lk;
        lk = lka + 2*N_train*t + N_train-1;
        mean = lk[root].mean;
        var = lk[root].var;
        peak = lk[root].peak; 
        infp = lk[root].infp;
      }
      v0 = (st[dt].d_SD[t]*st[dt].d_SD[t]) * st[dt].divt[-root];
      if (var>0) peak += 0.5*log(var);
      if (var+v0>0)
      { *logprob += peak - 0.5*log(var+v0) - 0.5*mean*mean/(var+v0);
      }
      else 
      { if (mean==0) 
        { *infprob += 1;
        }     
        else 
        { *logprob -= INFINITY;
        }
      }
      *infprob += infp;
    }
  }

  if (isnan(*logprob) || *logprob==INFINITY) abort();
}


/* COMPUTE LIKELIHOOD FUNCTION FOR PARENT OF NODE.  Finds the Gaussian 
   likelihood function for the location of the parent of a node in a tree, with
   respect to a particular target variable.  If the node is non-terminal,
   the likelihood depends on its location (which must be known).  If the node
   is terminal, the likelihood depends on the latent value associated with a
   it, or the data value, if there are no latent values (in which case the 
   data must be Gaussian).  For terminal nodes, if there is more than one tree,
   the likelihood also depends on the locations of the parents of this terminal
   node in the other trees, which must be known. 

   The likelihood function is Gaussian in shape, and is described by a "mean"
   and "variance", which are returned via the final two arguments. */

void dft_node_likelihood
( dft_spec *dft,		/* Diffusion tree model specification */
  int dt,			/* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  double *latent,		/* Latent vectors, or null */
  int x,			/* Index of node (+ve for terminals, -ve o.w.)*/
  int t,			/* Index of target variable (from 0) */
  double inv_temp,		/* Inverse temperature */
  double *mean,			/* Place to store "mean" of likelihood */
  double *var			/* Place to store "variance" of likelihood */
)
{
  int N_targets, px;
  double sd;

  N_targets = dft->N_targets;
  px = st[dt].parents[x];
  if (px>=0) abort();

  sd = st[dt].d_SD[t];

  if (x<0)
  { 
    if (st[dt].locations==0) abort();
    *mean = st[dt].locations[N_targets*(-x-1) + t];
    *var = (sd*sd) * (st[dt].divt[-x] - st[dt].divt[-px]);
  }
  else
  {
    dft_terminal_likelihood (dft, dt, st, latent, x, t, inv_temp, mean, var);
    *var += (sd*sd) * (1 - st[dt].divt[-px]);
  }

  if (*var<0) abort();
}


/* FIND LIKELIHOOD FUNCTION FOR A TERMINAL NODE.  Finds the Gaussian likelihood
   function for a terminal node, with respect to a particular target variable.
   The likelihood depends on the latent value associated with the terminal
   node, or on the data value, if there are no latent values (in which case
   the data must be Gaussian).  If there is more than one tree, the likelihood
   also depends on the locations of the parents of this terminal node in other
   trees, which must be known. 

   The "var" parameter can be set to zero, indicating an infinitely peaked
   likelihood at "mean". */

void dft_terminal_likelihood  
( dft_spec *dft,		/* Diffusion tree model specification */
  int dt,			/* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  double *latent,		/* Latent vectors, or null */
  int x,			/* Index of node (+ve for terminals, -ve o.w.)*/
  int t,			/* Index of target variable (from 0) */
  double inv_temp,		/* Inverse temperature */
  double *mean,			/* Place to store "mean" of likelihood */
  double *var			/* Place to store "variance" of likelihood */
)
{
  int N_targets;
  int dt2, px2;
  double sd2;

  N_targets = dft->N_targets;

  if (latent)
  { *mean = latent[N_targets*(x-1) + t];
    *var = 0;
  }
  else
  { *mean = train_targets[N_targets*(x-1) + t];
    *var = st->noise ? st->noise[t]*st->noise[t] / inv_temp : 0;
  }
  
  for (dt2 = 0; dt2<dft->N_trees; dt2++)
  { if (dt2!=dt)
    { if (st[dt2].locations==0) abort();
      px2 = st[dt2].parents[x];
      sd2 = st[dt2].d_SD[t];
      *var += (sd2*sd2) * (1 - st[dt2].divt[-px2]);
      *mean -= st[dt2].locations [(-px2-1)*N_targets + t];
    }
  }

  if (*var<0) abort();
}


/* COMPUTE DATA LIKELIHOOD FOR TREE.  Used only when node locations are not
   present.  Finds the Gaussian likelihood function pertaining to a particular 
   variable within the tree or subtree headed by the given root, returning it 
   via the mean, var, and peak parameters.  If the lk parameter is non-null, 
   the Gaussian likelihood functions for all the nodes in this sub-tree are 
   stored there, using indexes from -(N_train-1) to N_train (ie, the pointer 
   passed must be offset to allow this).  If the root parameter is zero, the 
   likelihood for the entire tree is found. 

   The "peak" is the log of the height of the likelihood function at "mean", 
   unless "var" is zero, in which case the -0.5*log(var) term is omitted.
   The values of "infp" is a count of points of infinite density.  A "peak" 
   value of +INFINITY should therefore not occur.  A "peak" value of -INFINITY 
   is  possible, from combining infinite peaks at different locations. */

void dft_tree_likelihood
( dft_spec *dft,		/* Diffusion tree model specification */
  int dt,			/* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  double *latent,		/* Latent vectors, or null */
  int t,			/* Variable to compute likelihood for */
  int root,			/* Root of sub-tree to find likelihood for */
  double inv_temp,		/* Inverse temperature */
  double *mean,			/* Place to store "mean" of likelihood */
  double *var,			/* Place to store "variance" of likelihood */
  double *peak,			/* Place to store log of peak value */
  int *infp,			/* Count of points of infinite density */
  dft_likelihood *lk		/* Place to store all likelihoods, or null */
)
{
  double m0, m1, var0, var1, peak0, peak1;
  int child0, child1, infp0, infp1;

  if (root<-(N_train-1) || root>N_train) abort();

  if (root==0)
  { for (root = 1; st[dt].parents[root]!=0; root = st[dt].parents[root]) ;
  }

  if (root>0)
  { dft_terminal_likelihood (dft, dt, st, latent, root, t, inv_temp, mean, var);
    if (*var>0)
    { *peak = - log_sqrt_2pi - 0.5*log(*var);
    }
    else
    { *peak = - log_sqrt_2pi;
    }
    *infp = 0;
  }
  else
  { child0 = chld(st[dt].nodes[-root],0);
    child1 = chld(st[dt].nodes[-root],1);
    if (child0==0 || child1==0) abort();
    dft_tree_likelihood (dft, dt, st, latent, t, child0, inv_temp,
                         &m0, &var0, &peak0, &infp0, lk);
    dft_tree_likelihood (dft, dt, st, latent, t, child1, inv_temp,
                         &m1, &var1, &peak1, &infp1, lk);
    dft_combine_likelihoods (child0, m0, var0, peak0, infp0, 
                             child1, m1, var1, peak1, infp1,
                             root, st[dt].divt, st[dt].d_SD[t], 
                             mean, var, peak, infp);
  }

  if (lk) 
  { lk[root].mean = *mean; 
    lk[root].var  = *var; 
    lk[root].peak = *peak;
    lk[root].infp = *infp;
  }
}


/* UPDATE TREE LIKELIHOOD.  Used only when node locations are not present.
   Updates the likelihoods for nodes recorded in the lk array (which has 
   likelihoods for each node for each variable) to account for changes 
   regarding non-terminal node x and/or its ancestors.  If the node passed
   is zero, nothing is done. */

void dft_update_likelihoods
( dft_spec *dft,		/* Diffusion tree model specification */
  int dt,			/* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  double *latent,		/* Latent vectors, or null */
  int x,			/* Node to update likelihood for (+ ancestors)*/
  dft_likelihood *lka		/* Likelihoods to update */
)
{
  double m0, m1, var0, var1, peak0, peak1, mean, var, peak;
  int child0, child1, infp0, infp1;
  double *divt, sd;
  dft_likelihood *lk;
  int t, y;

  if (x==0) return;

  if (x>0) abort();

  for (t = 0; t<dft->N_targets; t++)
  { 
    lk = lka + 2*N_train*t + N_train-1;

    sd = st[dt].d_SD[t];
    divt = st[dt].divt;

    for (y = x; y!=0; y = st[dt].parents[y])
    { child0 = chld(st[dt].nodes[-y],0);
      child1 = chld(st[dt].nodes[-y],1);
      m0    = lk[child0].mean;  var0  = lk[child0].var; 
      peak0 = lk[child0].peak;  infp0 = lk[child0].infp;
      m1    = lk[child1].mean;  var1  = lk[child1].var; 
      peak1 = lk[child1].peak;  infp1 = lk[child1].infp;
      dft_combine_likelihoods (child0, m0, var0, peak0, infp0, 
                               child1, m1, var1, peak1, infp1,
                               y, divt, sd, &lk[y].mean, &lk[y].var, 
                               &lk[y].peak, &lk[y].infp);
    }
  }
}


/* GET LIKELIHOOD FUNCTION FOR PARENT BY COMBINING THOSE OF ITS CHILDREN.

   Likelihoods are Gaussian-shaped, defined by "mean" and "var" values.  The
   "peak" values is the log of the height of the likelihood function at "mean", 
   unless "var" is zero, in which case the -0.5*log(var) term is omitted.  A
   "peak" value of +INFINITY should therefore not occur.  A "peak" value of 
   -INFINITY will result if two likelihoods with infinite peaks at different 
   locations are combined. */

void dft_combine_likelihoods
( int child0,		/* Variables pertaining to one child */
  double m0,
  double var0,
  double peak0,
  int infp0,
  int child1,		/* Variables pertaining to the other child */
  double m1,
  double var1,
  double peak1,
  int infp1,
  int parent,		/* Parent node */
  double *divt,		/* Array of divergence times */
  double sd,		/* Diffusion standard deviation */
  double *mean,		/* Variables to return combined likelihood in */
  double *var,
  double *peak,
  int *infp
)
{ 
  double ev0, ev1, var_ev0, var_ev1;

  ev0 = (sd*sd) * ((child0>0 ? 1 : divt[-child0]) - divt[-parent]);
  ev1 = (sd*sd) * ((child1>0 ? 1 : divt[-child1]) - divt[-parent]);
  if (ev0<0 || ev1<0) abort();

  var_ev0 = var0+ev0;
  var_ev1 = var1+ev1;

  if (var_ev0==INFINITY)
  { *mean = m1;
    *var = var1;
    *peak = peak1;
    *infp = infp1;
  }
  else if (var_ev1==INFINITY)
  { *mean = m0;
    *var = var0;
    *peak = peak0;
    *infp = infp0;
  }
  else if (var_ev0+var_ev1==0)
  { *mean = (m0+m1)/2;
    *var = 0;
    *peak = peak0 + peak1;
    *infp = infp0 + infp1;
    if (m0!=m1) *peak -= INFINITY;
    else *infp += 1;
  }
  else
  { *mean = var_ev0==0 ? m0 : var_ev1==0 ? m1
          : (m0*var_ev1 + m1*var_ev0) / (var_ev0+var_ev1);
    *var = var_ev0*var_ev1 / (var_ev0+var_ev1);
    if (var0>0) peak0 += 0.5*log(var0);
    if (var1>0) peak1 += 0.5*log(var1);
    *peak = peak0 + peak1;
    *infp = infp0 + infp1;
    if (var_ev0>0)
    { *peak += - 0.5*log(var_ev0) - 0.5*(*mean-m0)*(*mean-m0)/var_ev0;
    }
    if (var_ev1>0)
    { *peak += - 0.5*log(var_ev1) - 0.5*(*mean-m1)*(*mean-m1)/var_ev1;
    }
  }

  if (*var<0 || isnan(*mean) || isnan(*var) || isnan(*peak) || *peak==INFINITY)
  { abort();
  }
}
