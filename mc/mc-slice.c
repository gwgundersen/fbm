/* MC-SLICE.C - Procedures for performing slice sampling updates. */

/* Copyright (c) 1996 by Radford M. Neal 
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


/* LOCAL PROCEDURES. */

static void step_out   (mc_dynamic_state *, mc_iter *, int, double, double, int,
                        double *, double *);

static void pick_value (mc_dynamic_state *, mc_iter *, int, double, double, 
                        double, double);


/* PERFORM SLICE SAMPLING UPDATES FOR COORDINATES IN SOME RANGE. */

void mc_slice
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int firsti,		/* Index of first component to update (-1 for all) */
  int lasti,		/* Index of last component to update */
  int max_steps		/* Maximum number of intervals, zero for unlimited */
)
{
  double curr_q, slice_point, low_bnd, high_bnd;
  int k;

  if (firsti==-1) 
  { firsti = 0;
    lasti = ds->dim-1;
  }

  if (lasti>=ds->dim-1) lasti = ds->dim-1;
  if (firsti>lasti) firsti = lasti;
  
  for (k = firsti; k<=lasti; k++)
  {
    it->slice_calls += 1;

    if (!ds->know_pot)
    { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
      it->slice_evals += 1;
    }

    slice_point = ds->pot_energy + rand_exp();
    curr_q = ds->q[k];

    step_out (ds, it, k, slice_point, curr_q, max_steps, &low_bnd, &high_bnd);
    pick_value (ds, it, k, slice_point, curr_q, low_bnd, high_bnd);

    ds->know_pot = 1;
  }

  ds->know_grad = 0;
}


/* PERFORM OVERRELAXED SLICE SAMPLING UPDATES FOR COORDINATES IN SOME RANGE. */

void mc_slice_over
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int refinements,	/* Number of refinements to do to endpoints */
  float refresh_prob,	/* Probability of doing a refresh update */
  int firsti,		/* Index of first component to update (-1 for all) */
  int lasti,		/* Index of last component to update */
  int max_steps		/* Maximum number of intervals, zero for unlimited */
)
{
  double sf, curr_q, slice_point, low_bnd, high_bnd, olow_bnd, ohigh_bnd, width;
  int k, i, r;

  sf = it->stepsize_factor;

  if (firsti==-1) 
  { firsti = 0;
    lasti = ds->dim-1;
  }

  if (lasti>=ds->dim-1) lasti = ds->dim-1;
  if (firsti>lasti) firsti = lasti;
  
  for (k = firsti; k<=lasti; k++)
  {
    it->slice_calls += 1;

    if (!ds->know_pot)
    { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
      it->slice_evals += 1;
    }

    slice_point = ds->pot_energy + rand_exp();
    curr_q = ds->q[k];

    step_out (ds, it, k, slice_point, curr_q, max_steps, &low_bnd, &high_bnd);

    if (rand_uniform() < refresh_prob)
    { 
      pick_value (ds, it, k, slice_point, curr_q, low_bnd, high_bnd);
      ds->know_pot = 1;

      continue;
    }

    width = sf * ds->stepsize[k];
    r = refinements;

    if (high_bnd-low_bnd <= width*1.1)
    {
      while (r>0)
      {
        width /= 2;
        r -= 1;

        ds->q[k] = (low_bnd + high_bnd) / 2;

        mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
        it->slice_evals += 1;

        if (ds->pot_energy <= slice_point) break;

        if (ds->q[k]<curr_q) 
        { low_bnd = ds->q[k];
        }
        else 
        { high_bnd = ds->q[k];
        }
      }
    }

    olow_bnd = low_bnd;
    ohigh_bnd = high_bnd;

    while (r>0)
    { 
      width /= 2;

      ds->q[k] = low_bnd + width;

      mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
      it->slice_evals += 1;
      if (ds->pot_energy > slice_point)
      { low_bnd = ds->q[k];
      }

      ds->q[k] = high_bnd - width;

      mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
      it->slice_evals += 1;
      if (ds->pot_energy > slice_point)
      { high_bnd = ds->q[k];
      }

      r -= 1;
    }

    ds->q[k] = low_bnd + high_bnd - curr_q;

    if (ds->q[k]>=olow_bnd && ds->q[k]<=ohigh_bnd)
    { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
      it->slice_evals += 1;
      ds->know_pot = 1;
    }

    it->proposals += 1;

    if (ds->q[k]<olow_bnd || ds->q[k]>ohigh_bnd
     || ds->pot_energy > slice_point)
    { ds->q[k] = curr_q;
      it->rejects += 1;
      ds->know_pot = 0;
    }
  }

  ds->know_grad = 0;
}


/* FIND INTERVAL AROUND PART OF SLICE BY STEPPING OUT. */

static void step_out 
( mc_dynamic_state *ds,	/* Current state */
  mc_iter *it,		/* Description of this iteration */
  int k,		/* Index of coordinate being updated */
  double slice_point,	/* Potential energy level that defines slice */
  double curr_q,	/* Current value for coordinate */
  int max_steps,	/* Maximum number of intervals, zero for unlimited */
  double *low_bnd,	/* Place to store low bound for interval */
  double *high_bnd	/* Place to store high bound for interval */
)
{
  int low_steps, high_steps;
  int low_out, high_out;
  double sf;

  sf = it->stepsize_factor;

  low_steps  = max_steps==0 ? 1000000000 : rand_int(max_steps);
  high_steps = max_steps==0 ? 1000000000 : (max_steps-1) - low_steps;

  *low_bnd  = curr_q - rand_uniopen() * sf * ds->stepsize[k];
  *high_bnd = *low_bnd + sf * ds->stepsize[k];

  low_out = high_out = -1;

  for (;;)
  {
    if (low_out==-1)
    { ds->q[k] = *low_bnd;
      mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
      it->slice_evals += 1;
      low_out = ds->pot_energy > slice_point;
    }

    if (high_out==-1)
    { ds->q[k] = *high_bnd;
      mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
      it->slice_evals += 1;
      high_out = ds->pot_energy > slice_point;
    }

    if ((low_out || low_steps==0) && (high_out || high_steps==0)) return;

    if (!low_out && low_steps>0)
    { *low_bnd -= sf * ds->stepsize[k];
      low_out = -1;
    }
    
    if (!high_out && high_steps>0)
    { *high_bnd += sf * ds->stepsize[k];
      high_out = -1;
    }
  }
}


/* PICK VALUE FROM INTERVAL AROUND PART OF THE SLICE. */

static void pick_value
( mc_dynamic_state *ds,		/* Current state, updated with new value */
  mc_iter *it,			/* Description of this iteration */
  int k,			/* Index of coordinate being updated */
  double slice_point,		/* Potential energy level that defines slice */
  double curr_q,		/* Current value for coordinate */
  double low_bnd,		/* Low bound for interval */
  double high_bnd		/* High bound for interval */
)
{
  for (;;)
  {
    ds->q[k] = low_bnd + rand_uniopen() * (high_bnd-low_bnd);

    mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
    it->slice_evals += 1;

    if (ds->pot_energy<=slice_point) return;

    if (ds->q[k]<curr_q) 
    { low_bnd = ds->q[k];
    }
    else
    { high_bnd = ds->q[k];
    }
  }
}
