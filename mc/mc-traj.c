/* MC-TRAJ.C - Procedures for computing dynamical trajectories. */

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
#include "rand.h"
#include "log.h"
#include "mc.h"


/* HOW TO COMPUTE THE CURRENT TRAJECTORY. */

static mc_traj *tj;	/* Description of how to compute trajectory */
static mc_iter *it;	/* Description of iteration (contains approx order) */

static int na;		/* Number of distinct approximations used */
static int ta;		/* Total approximations, counting repeats for symmetry*/


/* SET UP FOR COMPUTING A TRAJECTORY.  Passed the description of how to
   compute the trajectory and the stepsizes to use.  Saves these, and
   sets other variables describing how the trajectory will be computed
   this time. */

void mc_traj_init
( mc_traj *tj0,		/* Description of how to compute trajectory */
  mc_iter *it0		/* Description of this iteration */
)
{ 
  tj = tj0;
  it = it0;

  if (tj->N_approx>0)
  { na = tj->N_approx;
    ta = na;
  }
  else
  { na = -tj->N_approx;
    ta = 2*na;
  }

  if (na>Max_approx)
  { fprintf (stderr, "Using too many energy approximations (max %d)\n",
             Max_approx);
    exit(1);
  }
}


/* PERMUTE ORDER OF APPROXIMATIONS. */

void mc_traj_permute (void)
{
  int i, j, t;

  for (i = 0; i<na; i++) 
  { it->approx_order[i] = i+1;
  }

  for (i = 0; i<na; i++)
  { t = it->approx_order[i];
    j = i + rand_int(na-i);
    it->approx_order[i] = it->approx_order[j];
    it->approx_order[j] = t;
  }
}


/* COMPUTE SOME NUMBER OF STEPS OF THE TRAJECTORY.  Updates the dynamical
   state by following the trajectory for some number of steps, computed
   as previously specified using mc_traj_init.  The number of steps may
   be negative, in which case the state is taken backwards along the
   trajectory that number of steps.  If the last parameter is non-zero,
   the potential energy will be evaluated at the last state, if it is
   convenient to do so. */

void mc_trajectory
( mc_dynamic_state *ds,	/* Dynamical state to update */
  int n,		/* Number of steps to compute */
  int need_pot		/* Need potential energy for last state? */
)
{
  double sf, sfh;
  int k, a, x, o;

  if (n==0) return;

  /* Set up stepsize, negating if a backward trajectory is desired.  Also
     divide stepsize by the number of approximations used, and multiply
     number of steps by the same number, to give a constant total time. */

  sf = it->stepsize_factor / ta;

  x = +1; 

  if (n<0)
  { n = -n;
    sf = -sf;
    x = -1;
  }

  n *= ta;

  sfh = sf / 2;

  /* Compute the trajectory, using as many approximations as asked for. */

  switch (tj->type)
  { 
    case 'L':
    { 
      if (tj->halfp)
      {
        /* Compute trajectory with initial and final half-steps for momentum .*/

        if (ds->know_grad!=it->approx_order[0])
        { mc_app_energy (ds, na, it->approx_order[0], 0, ds->grad);
        }
  
        for (k = 0; k<ds->dim; k++)
        { ds->p[k] -= sfh * ds->stepsize[k] * ds->grad[k];
        }

        a = 0;
  
        for (;;)
        {
          for (k = 0; k<ds->dim; k++)
          { ds->q[k] += sf * ds->stepsize[k] * ds->p[k];
          }
  
          a = (a+x+ta) % ta;
          o = a<na ? a : a==na ? 0 : 2*na-a;

          mc_app_energy (ds, na, it->approx_order[o], 
                         need_pot && n==1 ? &ds->pot_energy : 0, 
                         ds->grad);
  
          n -= 1;
          if (n==0) break;
  
          for (k = 0; k<ds->dim; k++)
          { ds->p[k] -= sf * ds->stepsize[k] * ds->grad[k];
          }
        }

        if (a!=0) abort();
  
        for (k = 0; k<ds->dim; k++)
        { ds->p[k] -= sfh * ds->stepsize[k] * ds->grad[k];
        }

        ds->know_grad = it->approx_order[0];
        ds->know_pot  = need_pot;
      }
      else 
      {
        /* Compute trajectory with initial and final half-steps for position. */

        for (k = 0; k<ds->dim; k++)
        { ds->q[k] += sfh * ds->stepsize[k] * ds->p[k];
        }

        a = x==-1 ? 0 : ta-1;
  
        for (;;)
        {
          a = (a+x+ta) % ta;
          o = a<na ? a : 2*na-a-1;

          mc_app_energy (ds, na, it->approx_order[o], 0, ds->grad);

          for (k = 0; k<ds->dim; k++)
          { ds->p[k] -= sf * ds->stepsize[k] * ds->grad[k];
          } 
  
          n -= 1;
          if (n==0) break;

          for (k = 0; k<ds->dim; k++)
          { ds->q[k] += sf * ds->stepsize[k] * ds->p[k];
          }
        }

        if (a != (x==-1 ? 0 : ta-1)) abort();
  
        for (k = 0; k<ds->dim; k++)
        { ds->q[k] += sfh * ds->stepsize[k] * ds->p[k];
        }

        ds->know_grad = 0;
        ds->know_pot  = 0;
      }
  
      break;
    }
  }

  /* Since momentum has changed, we don't know kinetic energy anymore. */

  ds->know_kinetic = 0;
}
