/* DFT-SAMPLE.C - Procedures that sample for latent vectors or node locations.*/

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
#include "prior.h"
#include "model.h"
#include "data.h"
#include "rand.h"
#include "ars.h"
#include "uars.h"
#include "phi.h"
#include "dft.h"
#include "dft-data.h"


/* SAMPLE FROM THE PRIOR OVER LOCATIONS AND LATENT VECTORS.  Ignores the
   data.  There must be only one tree. */

static double prior_sample
( dft_spec *dft,
  dft_state st,
  double *latent,
  int x,
  int t,
  double m,
  double v
)
{ 
  double ev0, ev1, value0, m1, v1, value1;
  int ch0, ch1;

  if (x>0)
  { 
    return latent[dft->N_targets*(x-1)+t] = m + sqrt(v)*rand_gaussian();
  } 
  else
  { 
    ch0 = chld(st->nodes[-x],0);
    ch1 = chld(st->nodes[-x],1);
    ev0 = (ch0>0 ? 1 : st->divt[-ch0]) - st->divt[-x];
    ev1 = (ch1>0 ? 1 : st->divt[-ch1]) - st->divt[-x];

    value0 = prior_sample (dft, st, latent, ch0, t, m, v+ev0);

    m1 = (m*ev0 + value0*v) / (v+ev0);
    v1 = v*ev0 / (v+ev0);
    value1 = prior_sample (dft, st, latent, ch1, t, m1, v1+ev1);

    m1 = (m1*ev1+value1*v1) / (v1+ev1);
    v1 = v1*ev1 / (v1+ev1);

    return st->locations[dft->N_targets*(-x-1)+t] = m1 + v1*rand_gaussian();
  }
}

void dft_sample_prior
( dft_spec *dft,		/* Specification for diffusion tree model */
  dft_state st,			/* Pointers to state, including locations */
  double *latent		/* Place to store sampled latent vectors */
)
{ 
  int root, t;

  if (dft->N_trees!=1) abort();

  for (root = 1; st->parents[root]!=0; root = st->parents[root]) ;

  for (t = 0; t<dft->N_targets; t++)
  {
    prior_sample (dft, st, latent, root, t, 0, root>0 ? 1 : st->divt[-root]);
  }
}


/* SAMPLE FOR LATENT VECTORS GIVEN NODE LOCATIONS.  Draws latent vectors 
   randomly from their conditional distribution given the data vectors and 
   the locations of non-terminal nodes. */

struct ars_info { int target; double pmean, pvar, inv_temp; };

static double ars_logp (double x, double *deriv, void *in)
{ struct ars_info *info = in;
  double l, d;
  l = - 0.5 * (x - info->pmean) * (x - info->pmean) / info->pvar;
  d = - (x - info->pmean) / info->pvar;
  if (info->target==1)
  { l += - info->inv_temp * log(1+exp(-x));
    d += info->inv_temp / (1+exp(x));
  }
  else if (info->target==0)
  { l += - info->inv_temp * log(1+exp(x));
    d += - info->inv_temp / (1+exp(-x));
  }
  else
  { abort();
  }
  *deriv = d;
  return l;
}

struct uars_info { double trans, scale, alpha, inv_temp; };

static double uars_logp (double y, void *in)
{ struct uars_info *info = in;
  double x;
  x = info->scale * Phi_inverse(y) + info->trans;
  return - info->inv_temp * ((info->alpha+1)/2) * log(1+x*x/info->alpha);
}

static double latent_sample
( dft_spec *dft,
  model_specification *m,
  dft_state st,
  double pmean,
  double pvar,
  double target,
  double inv_temp,
  int t
)
{
  double mean, var, nv;

  if (m->type=='R' && m->noise.alpha[2]==0)
  { if (inv_temp==0)
    { mean = pmean;
      var = pvar;
    }
    else
    { nv = st->noise[t] * st->noise[t] / inv_temp;
      mean = (target*pvar + pmean*nv) / (pvar + nv);
      var = pvar * nv / (pvar + nv);
    }
/*fprintf(stderr,
"latent_sample: inv_temp=%.3f pmean=%.3f pvar=%.3f target=%.3f noise=%.3f : %.3f %.3f\n",inv_temp,pmean,pvar,target,st->noise[t]*st->noise[t],mean,var);*/
    return mean + sqrt(var) * rand_gaussian();
  }

  else if (m->type=='R' && m->noise.alpha[2]!=0)
  { struct uars_info info;
    info.trans = (pmean-target) / st->noise[t];
    info.scale = sqrt(pvar) / st->noise[t];
    info.alpha = m->noise.alpha[2];
    info.inv_temp = inv_temp;
    return pmean + sqrt(pvar) * 
           Phi_inverse (uars(uars_logp, Phi((target-pmean)/sqrt(pvar)), &info));
  } 

  else if (m->type=='B')
  { struct ars_info info;
    info.target = target;
    info.pmean = pmean;
    info.pvar = pvar;
    info.inv_temp = inv_temp;
    return ars (0.0, 1.0, ars_logp, &info);
  }

  else /* unknown model type */
  { abort();
  }

}

void dft_sample_latent
( dft_spec *dft,		/* Specification for diffusion tree model */
  model_specification *m,	/* Specification for data model */
  dft_state st,			/* Pointers to state */
  double inv_temp,		/* Inverse temperature */
  double *latent		/* Place to store sampled latent vectors */
)
{ 
  int N_targets, N_trees, dt;

  if (!st->locations) abort();

  N_targets = dft->N_targets;
  N_trees = dft->N_trees;

  if (m==0)
  { int i, t;
    for (i = 0; i<N_train; i++)
    { for (t = 0; t<N_targets; t++)
      { latent[i*N_targets+t] = train_targets[i*N_targets+t];
      }
    }
  }

  else 
  { 
    double pmean, pvar;
    int t, x, px;

    for (t = 0; t<N_targets; t++)
    {
      for (x = 1; x<=N_train; x++)
      { 
        pmean = 0;
        pvar = 0;
        for (dt = 0; dt<N_trees; dt++)
        { px = st[dt].parents[x];
          pvar += (st[dt].d_SD[t] * st[dt].d_SD[t]) 
                    * (1 - (px==0 ? 0 : st[dt].divt[-px]));
          pmean += px==0 ? 0 : st[dt].locations[N_targets*(-px-1)+t];
        }

        latent[N_targets*(x-1)+t] = latent_sample (dft, m, st, pmean, pvar,
                                 train_targets[N_targets*(x-1)+t], inv_temp, t);
      }
    }
  }
}


/* GIBBS SAMPLING FOR LATENT VECTORS.  Does a Gibbs sampling scan over
   latent vectors, sampling from their conditional distributions given
   other latent vectors and the data (node locations are assumed absent). 
   There must be only one tree for this procedure to be used. 

   Uses the latent_sample procedure defined above for dft_sample_latent. */

static void latent_gibbs
( dft_spec *dft,
  model_specification *m,
  dft_state st,
  double *latent,
  int x,
  int t,
  double pmean,
  double pvar,
  double inv_temp,
  dft_likelihood *lk
)
{
  double omean, ovar, sd, d;
  int N_targets, c, ch[2];

  if (x<-(N_train-1) || x>N_train || x==0) abort();

  N_targets = dft->N_targets;

  if (x>0) /* Terminal node */
  {
    latent[N_targets*(x-1)+t] = latent_sample (dft, m, st, pmean, pvar,
                                 train_targets[N_targets*(x-1)+t], inv_temp, t);
    lk[x].mean = latent[N_targets*(x-1)+t];
    lk[x].var = 0;
    lk[x].peak = -log_sqrt_2pi;
  }

  else /* Non-terminal node */
  { 
    sd = st->d_SD[t];

    for (c = 0; c<=1; c++)
    { ch[c] = chld(st->nodes[-x],c);
    }

    d = st->divt[-x];

    for (c = 0; c<=1; c++)
    { omean = lk[ch[c]].mean;
      ovar  = lk[ch[c]].var + (sd*sd) * (ch[c]>0 ? 1-d : st->divt[-ch[c]]-d);
      omean = (omean*pvar + pmean*ovar) / (pvar+ovar);
      ovar  = pvar*ovar / (pvar+ovar);
      ovar += (sd*sd) * (ch[1-c]>0 ? 1-d : st->divt[-ch[1-c]]-d);
      latent_gibbs (dft, m, st, latent, ch[1-c], t, omean, ovar, inv_temp, lk);
    }

    dft_combine_likelihoods 
       (ch[0], lk[ch[0]].mean, lk[ch[0]].var, lk[ch[0]].peak, lk[ch[0]].infp,
        ch[1], lk[ch[1]].mean, lk[ch[1]].var, lk[ch[1]].peak, lk[ch[1]].infp,
        x, st->divt, sd,
        &lk[x].mean, &lk[x].var, &lk[x].peak, &lk[x].infp);
  }
} 

void dft_gibbs_latent
( dft_spec *dft,		/* Specification for diffusion tree model */
  model_specification *m,	/* Specification for data model */
  dft_state st,			/* Pointers to state */
  double inv_temp,		/* Inverse temperature */
  double *latent		/* The latent vectors */
)
{
  double mean, var, peak;
  int N_targets, N_trees, infp;
  dft_likelihood *lk0, *lk;
  int root, t;

  if (dft->N_trees>1) abort();

  if (N_train<2) return;

  N_targets = dft->N_targets;
  N_trees = dft->N_trees;

  if (m==0)
  { int i, t;
    for (i = 0; i<N_train; i++)
    { for (t = 0; t<N_targets; t++)
      { latent[i*N_targets+t] = train_targets[i*N_targets+t];
      }
    }
  }

  else 
  { 
    lk0 = chk_alloc (2*N_train, sizeof *lk0);
    lk = lk0 + N_train-1;

    for (root = 1; st->parents[root]!=0; root = st->parents[root]) ;

    for (t = 0; t<N_targets; t++)
    { 
      dft_tree_likelihood (dft, 0, st, latent, t, root, inv_temp,
                           &mean, &var, &peak, &infp, lk);

      latent_gibbs (dft, m, st, latent, root, t,
                    0, st->d_SD[t]*st->d_SD[t] * st->divt[-root], inv_temp, lk);
    }

    free(lk0);
  }
}


/* GIBBS SAMPLING FOR LOCATIONS OF NON-TERMINAL NODES.  Generates new locations
   for non-terminal nodes in all trees, one tree at at time, conditional on the
   current locations for the non-terminal nodes in the other trees, and the 
   latent vectors, if they are present, or the data vectors, if latent vectors 
   are not present.  Latent vectors must be present if the data is not real, or 
   if the data is real but with t-distributed noise. */

static void sample_loc 
( dft_spec *dft,
  dft_state st,
  int dt,
  int t,
  dft_likelihood *lk,
  int x,
  double m,
  double v
)
{ 
  double *loc;
  double vc[2], vp[2], mp[2];
  int ch[2];
  double mx, vx, sd;
  int c;
 
  loc = st[dt].locations + dft->N_targets*(-x-1) + t;
  sd = st[dt].d_SD[t];

  for (c = 0; c<=1; c++)
  { ch[c] = chld(st[dt].nodes[-x],c);
    vc[c] = (sd*sd) * ((ch[c]>0 ? 1 : st[dt].divt[-ch[c]]) - st[dt].divt[-x]);
    vp[c] = vc[c] + lk[ch[c]].var;
    mp[c] = lk[ch[c]].mean;
  }

  mx = (m*vp[0]*vp[1] + mp[0]*vp[1]*v + mp[1]*vp[0]*v) 
         / (vp[0]*vp[1] + vp[1]*v + vp[0]*v);

  vx = vp[0]*vp[1]*v / (vp[0]*vp[1] + vp[1]*v + vp[0]*v);

  *loc = mx + sqrt(vx)*rand_gaussian();
  if (isnan(*loc)) abort();

  for (c = 0; c<=1; c++)
  { if (ch[c]<0)
    { sample_loc (dft, st, dt, t, lk, ch[c], *loc, vc[c]);
    }
  }
}

void dft_gibbs_locations
( dft_spec *dft,		/* Specification for diffusion tree model */
  model_specification *m,	/* Specification for data model */
  dft_state st,			/* Pointers to state, including locations */
  double inv_temp,		/* Inverse temperature */
  double *latent		/* Latent vectors for training cases, or null */
)
{ 
  dft_likelihood *lk0, *lk;
  double mean, var, peak;
  int root, dt, infp, t;

  if (latent==0 && m!=0 && (m->type!='R' || m->noise.alpha[2]!=0))
  { fprintf(stderr,
"Can't sample for node locations with no latent vectors if model is non-Gaussian\n");
    exit(1);
  }

  if (N_train<=1) return;
 
  lk0 = chk_alloc (2*N_train, sizeof *lk);
  lk = lk0 + N_train-1;

  for (dt = 0; dt<dft->N_trees; dt++)
  { for (root = 1; st[dt].parents[root]!=0; root = st[dt].parents[root]) ;
    for (t = 0; t<dft->N_targets; t++)
    { dft_tree_likelihood (dft, dt, st, latent, t, root, inv_temp,
                           &mean, &var, &peak, &infp, lk);
      sample_loc (dft, st, dt, t, lk, root, 0, 
                  st[dt].divt[-root] * st[dt].d_SD[t]*st[dt].d_SD[t]);
    }
  }

  free(lk0);
}


/* SAMPLE FOR LOCATIONS OF TERMINAL NODES.  Generate locations for terminal 
   nodes, given latent vectors and locations for the parents of the terminal 
   nodes.  If there is only one tree, the terminal node location is just equal 
   to the latent value.  If there is more than one tree, the terminal locations
   must be randomly drawn from their distribution constrained to the hyperplane
   on which they add up to the latent value. */

void dft_sample_terminals
( dft_spec *dft,		/* Specification for diffusion tree model */
  dft_state st,			/* Pointers to state */
  double *latent,		/* Latent vectors for training cases */
  double *terminals[Max_trees]	/* Places to store terminal locations */
)
{
  double m[Max_trees], s[Max_trees], w[Max_trees], u[Max_trees], r[Max_trees];
  int N_targets, N_trees;
  double v, ul, ur, uw;
  int dt, x, px, t, c;

  N_targets = dft->N_targets;
  N_trees = dft->N_trees;

  if (N_trees==1)
  { 
    for (x = 1; x<=N_train; x++)
    { for (t = 0; t<N_targets; t++)
      { terminals[0][N_targets*(x-1)+t] = latent[N_targets*(x-1)+t];
      }
    }
  }

  else /* more than one tree */
  { 
    for (x = 1; x<=N_train; x++)
    { for (t = 0; t<N_targets; t++)
      { 
        /* Find the shift and scale transformations needed to convert the
           Gaussian prior for the terminal nodes of the trees based on their
           parents to a spherical Gaussian centred at the origin. */

        for (dt = 0; dt<N_trees; dt++)
        { px = st[dt].parents[x];
          m[dt] = px==0 ? 0 : st[dt].locations[N_targets*(-px-1)+t];
          s[dt] = st[dt].d_SD[t] * (px==0 ? 1 : sqrt(1-st[dt].divt[-px]));
        }

        /* Transform a point that lies on the constraint plane (where the 
           sum of the terminals is the actual latent value), taking care
           to avoid a point that would transform to infinity. */

        c = 0;
        for (dt = 0; dt<N_trees; dt++)
        { if (s[dt]!=0) c += 1;
        }

        if (c==0) abort();
        
        v = latent[N_targets*(x-1)+t] / c;

        for (dt = 0; dt<N_trees; dt++)
        { w[dt] = s[dt]==0 ? 0 : (v - m[dt]) / s[dt];
        }

        /* Find a unit vector normal to the transformed constraint plane. */

        ul = 0;
        for (dt = 0; dt<N_trees; dt++)
        { ul += s[dt]*s[dt];
        }
        ul = sqrt(ul);

        for (dt = 0; dt<N_trees; dt++)
        { u[dt] = s[dt]/ul;
        }

        /* Generate a random vector from the spherical Gaussian. */

        for (dt = 0; dt<N_trees; dt++)
        { r[dt] = rand_gaussian();
        }

        /* Project the random vector onto the constraint plane. */

        ur = uw = 0;
        for (dt = 0; dt<N_trees; dt++)
        { ur += u[dt]*r[dt];
          uw += u[dt]*w[dt];
        }

        for (dt = 0; dt<N_trees; dt++)
        { r[dt] += (uw-ur)*u[dt];
        }

        /* Set the terminals to the result of transforming this value
           back to the original coordinate system. */

        for (dt = 0; dt<N_trees; dt++)
        { terminals[dt][N_targets*(x-1)+t] = m[dt] + s[dt]*r[dt];
          if (isnan(terminals[dt][N_targets*(x-1)+t])) abort();
        }
      }
    }
  }
}
