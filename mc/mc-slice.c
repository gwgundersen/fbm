/* MC-SLICE.C - Procedures for performing slice sampling updates. */

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


/* LOCAL PROCEDURES. */

static void step_out   (mc_dynamic_state *, mc_iter *, int, double, double, int,
                        double *, double *);

static void dbl_out    (mc_dynamic_state *, mc_iter *, int, double, double, int,
                        double *, double *);

static int dbl_ok      (mc_dynamic_state *, mc_iter *, int, double, double,
                        double, double);

static void pick_value (mc_dynamic_state *, mc_iter *, int, double, double, int,
                        double, double, int, double);


/* PERFORM ONE-VARIABLE SLICE SAMPLING UPDATES FOR COORDINATES IN SOME RANGE. */

void mc_slice_1
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int firsti,		/* Index of first component to update (-1 for all) */
  int lasti,		/* Index of last component to update */
  int max_steps,	/* Maximum number of intervals, zero for unlimited;
			   if negative, intervals are found by doubling */
  int r_update,		/* Update just one component at random? */
  int s_factor,		/* Factor for faster shrinkage */
  double s_threshold	/* Threshold for faster shrinkage */
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

  if (r_update)
  { firsti = lasti = firsti + (int)(rand_uniform()*(lasti-firsti+1));
  }
  
  for (k = firsti; k<=lasti; k++)
  {
    it->slice_calls += 1;

    if (!ds->know_pot)
    { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
      it->slice_evals += 1;
    }

    slice_point = ds->pot_energy + rand_exp();
    curr_q = ds->q[k];

    if (max_steps>=0)
    { step_out (ds, it, k, slice_point, curr_q, max_steps, &low_bnd, &high_bnd);
    }
    else
    { dbl_out (ds, it, k, slice_point, curr_q, -max_steps, &low_bnd, &high_bnd);
    }

    pick_value (ds, it, k, slice_point, curr_q, max_steps, low_bnd, high_bnd,
                s_factor, s_threshold);
  }

  ds->know_grad = 0;
}


/* PERFORM MULTIVARIATE SLICE SAMPLING WITH HYPERRECTANGLES. */

void mc_slice
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  mc_value *save,	/* Place to save current state */
  mc_value *lowb,	/* Storage for low bounds of hyperrectangle */
  mc_value *highb,	/* Storage for high bounds of hyperrectangle */
  int g_shrink		/* Shrink based on gradient? */
)
{
  double init_energy, slice_point, sf, maxp, pr;
  int k, maxk;

  it->slice_calls += 1;

  if (!ds->know_pot)
  { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
    ds->know_pot = 1;
    it->slice_evals += 1;
  }

  slice_point = ds->pot_energy + rand_exp();
  init_energy = ds->pot_energy;

  mc_value_copy (save, ds->q, ds->dim);

  sf = it->stepsize_factor;

  for (k = 0; k<ds->dim; k++) 
  { lowb[k]  = ds->q[k] - sf * ds->stepsize[k] * rand_uniopen();
    highb[k] = lowb[k] + sf * ds->stepsize[k];
  }
 
  for (;;)
  { 
    for (k = 0; k<ds->dim; k++) 
    { ds->q[k]  = lowb[k] + rand_uniopen() * (highb[k]-lowb[k]);
    }

    mc_app_energy (ds, 1, 1, &ds->pot_energy, g_shrink ? ds->grad : 0);
    ds->know_pot = 1;
    ds->know_grad = g_shrink;
    it->slice_evals += 1;

    if (ds->pot_energy<=slice_point)
    { it->delta = ds->pot_energy - init_energy;
      return;
    }

    if (g_shrink!=0)
    { 
      maxp = -1;

      for (k = 0; k<ds->dim; k++) 
      { 
        pr = ds->grad[k] * (highb[k]-lowb[k]);
        if (pr<0) pr = -pr;

        if (pr>maxp)
        { maxp = pr;
          maxk = k;
        }
      }
    }

    for (k = 0; k<ds->dim; k++) 
    { 
      if (g_shrink==2) 
      { pr = ds->grad[k] * (highb[k]-lowb[k]);
        if (pr<0) pr = -pr;
      }

      if (g_shrink==0 || k==maxk || g_shrink==2 && pr>=maxp/2)
      { if (ds->q[k]>save[k])
        { highb[k] = ds->q[k];
        }
        else
        { lowb[k] = ds->q[k];
        }
      }
    }
  }
}


/* PERFORM MULTIVARIATE GAUSSIAN SLICE SAMPLING. */

void mc_slice_gaussian
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  mc_value *save,	/* Place to save current state */
  mc_value *wsum,	/* Place to store weighted sum of crumbs */
  int e_shrink		/* Shrink dist. based on energy of trial points? */
)
{
  double init_energy, slice_point, sf, w, wo, totw;
  int k;

  it->slice_calls += 1;

  if (!ds->know_pot)
  { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
    ds->know_pot = 1;
    it->slice_evals += 1;
  }

  slice_point = ds->pot_energy + rand_exp();  
  init_energy = ds->pot_energy;

  mc_value_copy (save, ds->q, ds->dim);

  sf = it->stepsize_factor;
  w = 1/(sf*sf);

  for (k = 0; k<ds->dim; k++)
  { wsum[k] = w * (save[k] + rand_gaussian() * ds->stepsize[k] / sqrt(w));
  }

  totw = w;

  for (;;)
  { 
    for (k = 0; k<ds->dim; k++) 
    { ds->q[k]  = wsum[k]/totw + rand_gaussian() * ds->stepsize[k] / sqrt(totw);
    }

    mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
    ds->know_pot = 1;
    ds->know_grad = 0;
    it->slice_evals += 1;

    if (ds->pot_energy<=slice_point)
    { it->delta = ds->pot_energy - init_energy;
      return;
    }
  
    if (e_shrink)
    { wo = w;
      w = totw * (ds->pot_energy - slice_point - 1);
      if (w<wo) w = wo;
    }

    for (k = 0; k<ds->dim; k++) 
    { wsum[k] += w * (save[k] + rand_gaussian() * ds->stepsize[k] / sqrt(w));
    }

    totw += w;
  }
}


/* PERFORM OVERRELAXED SLICE SAMPLING UPDATES FOR COORDINATES IN SOME RANGE. */

void mc_slice_over
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int refinements,	/* Number of refinements to do to endpoints */
  float refresh_prob,	/* Probability of doing a refresh update */
  int firsti,		/* Index of first component to update (-1 for all) */
  int lasti,		/* Index of last component to update */
  int max_steps,	/* Maximum number of intervals, zero for unlimited */
  int r_update		/* Update just one component at random? */
)
{
  double low_bnd, high_bnd, olow_bnd, ohigh_bnd, dlow_bnd, dhigh_bnd;
  double sf, curr_q, slice_point, width;
  int k, r;

  sf = it->stepsize_factor;

  if (firsti==-1) 
  { firsti = 0;
    lasti = ds->dim-1;
  }

  if (lasti>=ds->dim-1) lasti = ds->dim-1;
  if (firsti>lasti) firsti = lasti;

  if (r_update)
  { firsti = lasti = firsti + (int)(rand_uniform()*(lasti-firsti+1));
  }
  
  for (k = firsti; k<=lasti; k++)
  {
    it->slice_calls += 1;

    if (!ds->know_pot)
    { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
      it->slice_evals += 1;
    }

    slice_point = ds->pot_energy + rand_exp();
    curr_q = ds->q[k];

    if (max_steps>=0)
    { step_out (ds, it, k, slice_point, curr_q, max_steps, &low_bnd, &high_bnd);
    }
    else
    { dbl_out (ds, it, k, slice_point, curr_q, -max_steps, &low_bnd, &high_bnd);
    }

    if (rand_uniform() < refresh_prob)
    { pick_value (ds, it, k, slice_point, curr_q, max_steps, low_bnd, high_bnd,
                  0, 0.0);
      continue;
    }

    width = sf * ds->stepsize[k];
    r = refinements;

    if (max_steps<0) 
    { width = high_bnd-low_bnd;  /* So we'll shrink 'til midpoint inside */
      dlow_bnd = low_bnd;
      dhigh_bnd = high_bnd;
    }

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

    if (ds->q[k]<olow_bnd || ds->q[k]>ohigh_bnd || ds->pot_energy > slice_point
     || (max_steps<0 && !dbl_ok(ds,it,k,slice_point,curr_q,dlow_bnd,dhigh_bnd)))
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
  double width;

  width = it->stepsize_factor * ds->stepsize[k];

  low_steps  = max_steps==0 ? 1000000000 
             : max_steps==1 ? 0
             : rand_int(max_steps);
  high_steps = max_steps==0 ? 1000000000 
             : (max_steps-1) - low_steps;

  *low_bnd  = curr_q - rand_uniopen() * width;
  *high_bnd = *low_bnd + width;

  while (low_steps>0)
  {
    ds->q[k] = *low_bnd;
    mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
    it->slice_evals += 1;

    if (ds->pot_energy>slice_point) break;

    *low_bnd -= width;
    low_steps -= 1;
  }

  while (high_steps>0)
  {
    ds->q[k] = *high_bnd;
    mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
    it->slice_evals += 1;

    if (ds->pot_energy>slice_point) break;

    *high_bnd += width;
    high_steps -= 1;
  }
}


/* FIND INTERVAL AROUND PART OF SLICE BY DOUBLING. */

static void dbl_out
( mc_dynamic_state *ds,	/* Current state */
  mc_iter *it,		/* Description of this iteration */
  int k,		/* Index of coordinate being updated */
  double slice_point,	/* Potential energy level that defines slice */
  double curr_q,	/* Current value for coordinate */
  int max_int,		/* Maximum number of intervals to go through */
  double *low_bnd,	/* Place to store low bound for interval */
  double *high_bnd	/* Place to store high bound for interval */
)
{
  int low_out, high_out;
  double width;

  width = it->stepsize_factor * ds->stepsize[k];

  *low_bnd  = curr_q - rand_uniopen() * width;
  *high_bnd = *low_bnd + width;

  ds->q[k] = *low_bnd;
  mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
  it->slice_evals += 1;
  low_out = ds->pot_energy>slice_point;

  ds->q[k] = *high_bnd;
  mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
  it->slice_evals += 1;
  high_out = ds->pot_energy>slice_point;

  while (max_int>1 && (!low_out || !high_out))
  { 
    if (rand_int(2))
    { *low_bnd -= (*high_bnd - *low_bnd);
      ds->q[k] = *low_bnd;
      mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
      it->slice_evals += 1;
      low_out = ds->pot_energy>slice_point;
    }
    else
    { *high_bnd += (*high_bnd - *low_bnd);
      ds->q[k] = *high_bnd;
      mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
      it->slice_evals += 1;
      high_out = ds->pot_energy>slice_point;
    }

    max_int -= 1;
  }
}


/* CHECK THAT DOUBLING WOULD FIND SAME INTERVAL FROM NEW POINT.  Returns 1
   if it would, 0 if it would not.  Will set ds->know_pot to 0 if it 
   disturbs ds->pot_energy. */

static int dbl_ok
( mc_dynamic_state *ds,	/* Current state */
  mc_iter *it,		/* Description of this iteration */
  int k,		/* Index of coordinate being updated */
  double slice_point,	/* Potential energy level that defines slice */
  double curr_q,	/* Current value for coordinate */
  double low_bnd,	/* Low bound for interval */
  double high_bnd	/* High bound for interval */
)
{
  double width, mid;
  double new_q;
  int diff;

  width = it->stepsize_factor * ds->stepsize[k];

  new_q = ds->q[k];
  diff = 0;

  while (high_bnd-low_bnd>1.1*width)
  { mid = (high_bnd+low_bnd)/2;
    if (!diff) 
    { diff = (curr_q<mid) != (new_q<mid);
    }
    if (diff)
    { ds->q[k] = mid;
      mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
      it->slice_evals += 1;
      ds->know_pot = 0;
    }
    if (new_q<mid)
    { high_bnd = mid;
      ds->q[k] = low_bnd;
    }
    else
    { low_bnd = mid;
      ds->q[k] = high_bnd;
    }
    if (diff && ds->pot_energy>slice_point)
    { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
      it->slice_evals += 1;
      if (ds->pot_energy>slice_point)
      { ds->q[k] = new_q;
        return 0;
      }
    }
  }

  ds->q[k] = new_q;
  return 1;
}


/* PICK VALUE FROM INTERVAL AROUND PART OF THE SLICE.  Sets ds->know_pot
   appropriately. */

static void pick_value
( mc_dynamic_state *ds,		/* Current state, updated with new value */
  mc_iter *it,			/* Description of this iteration */
  int k,			/* Index of coordinate being updated */
  double slice_point,		/* Potential energy level that defines slice */
  double curr_q,		/* Current value for coordinate */
  int max_steps,		/* Negative if interval was found by doubling */
  double low_bnd,		/* Low bound for interval */
  double high_bnd,		/* High bound for interval */
  int s_factor,			/* Factor for faster shrinkage */
  double s_threshold		/* Threshold for faster shrinkage */
)
{
  double nlow_bnd, nhigh_bnd, range;

  nlow_bnd = low_bnd;
  nhigh_bnd = high_bnd;

  for (;;)
  {
    ds->q[k] = nlow_bnd + rand_uniopen() * (nhigh_bnd-nlow_bnd);

    mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);    
    it->slice_evals += 1;
    ds->know_pot = 1;

    if (ds->pot_energy<=slice_point && (max_steps>=0 
         || dbl_ok(ds,it,k,slice_point,curr_q,low_bnd,high_bnd)))
    { return;
    }

    if (s_factor>=0)
    { 
      if (ds->q[k]<curr_q) 
      { nlow_bnd = ds->q[k];
      }
      else
      { nhigh_bnd = ds->q[k];
      }
    }

    if ((s_factor<-1 || s_factor>1) && ds->pot_energy>slice_point+s_threshold)
    { 
      range = (nhigh_bnd-nlow_bnd) / (s_factor>0 ? s_factor : -s_factor);
      if (range<=0) abort();
      while (nlow_bnd+range<curr_q) 
      { nlow_bnd += range;
      }
      while (nhigh_bnd-range>curr_q) 
      { nhigh_bnd -= range;
      }
    } 
  }
}


/* PERFORM MULTIVARIATE SLICE SAMPLING WITH INSIDE REFLECTION. */

void mc_slice_inside
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int steps,		/* Number of steps to take */
  mc_value *q_save,	/* Place to save q values */
  mc_value *p_save	/* Place to save p values (used here for grad) */
)
{
  double slice_point;  
  double gmag, proj;
  double old_pot;
  double sf;

  int rejects, rejected;
  int k, j;

  sf = it->stepsize_factor;

  if (!ds->know_pot || !ds->know_grad)
  { mc_app_energy (ds, 1, 1, &ds->pot_energy, ds->grad);
    ds->know_pot = 1;
    ds->know_grad = 1;
  }

  slice_point = ds->pot_energy + rand_exp();

  rejects = 0;
  
  for (k = 0; k<steps; k++)
  {
    mc_value_copy (q_save, ds->q, ds->dim);
    mc_value_copy (p_save, ds->grad, ds->dim);
    old_pot = ds->pot_energy;

    for (j = 0; j<ds->dim; j++) 
    { ds->q[j] += sf * ds->stepsize[j] * ds->p[j];
    }

    mc_app_energy (ds, 1, 1, &ds->pot_energy, ds->grad);

    if (ds->pot_energy>slice_point) 
    { 
      mc_value_copy (ds->q, q_save, ds->dim);
      mc_value_copy (ds->grad, p_save, ds->dim);
      ds->pot_energy = old_pot;

      gmag = 0;
      for (j = 0; j<ds->dim; j++) 
      { gmag += ds->grad[j]*ds->grad[j];
      }

      proj = 0;
      for (j = 0; j<ds->dim; j++) 
      { proj += ds->p[j]*ds->grad[j];
      }

      rejected = 1;

      if (gmag>1e-30)
      { 
        for (j = 0; j<ds->dim; j++)
        { ds->q[j] += 
            sf * ds->stepsize[j] * (2*ds->grad[j]*proj/gmag - ds->p[j]);
        }      

        mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);

        if (ds->pot_energy>slice_point)
        { for (j = 0; j<ds->dim; j++)
          { ds->p[j] = 2*ds->grad[j]*proj/gmag - ds->p[j];
          }      
          rejected = 0;
        }

        mc_value_copy (ds->q, q_save, ds->dim);
        ds->pot_energy = old_pot;
      }

      for (j = 0; j<ds->dim; j++) ds->p[j] = - ds->p[j];
     
      rejects += rejected;
    }
  }

  for (j = 0; j<ds->dim; j++) ds->p[j] = - ds->p[j];

  it->proposals += k;
  it->rejects += rejects;
}


/* PERFORM MULTIVARIATE SLICE SAMPLING WITH OUTSIDE REFLECTION. */

void mc_slice_outside
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int steps,		/* Number of steps to take */
  int in_steps,		/* Max number of steps inside slice */
  mc_value *q_save,	/* Place to save q values */
  mc_value *p_save	/* Place to save p values (used here for grad) */
)
{
  double slice_point;  
  double gmag, proj;
  double old_pot;
  double sf;

  int j, k, ki;

  sf = it->stepsize_factor;

  if (!ds->know_pot || !ds->know_grad)
  { mc_app_energy (ds, 1, 1, &ds->pot_energy, ds->grad);
    ds->know_pot = 1;
    ds->know_grad = 1;
  }

  mc_value_copy (q_save, ds->q, ds->dim);
  mc_value_copy (p_save, ds->grad, ds->dim);
  old_pot = ds->pot_energy;

  slice_point = ds->pot_energy + rand_exp();

  /*fprintf(stderr,"E: %f %f",ds->pot_energy,slice_point);*/
  /*printf("\n%f %f\n",ds->q[0],ds->q[1]);*/

  k = 0;
  ki = 0;

  for (;;)
  {
    for (j = 0; j<ds->dim; j++) 
    { ds->q[j] += sf * ds->stepsize[j] * ds->p[j];
    }

    /*printf("%f %f\n",ds->q[0],ds->q[1]);*/

    k += 1;

    mc_app_energy (ds, 1, 1, &ds->pot_energy, ds->grad);

    if (ds->pot_energy<=slice_point) ki += 1;

    if (k>=steps || ki>=in_steps) break;

    if (ds->pot_energy>slice_point) 
    { 
      gmag = 0;
      for (j = 0; j<ds->dim; j++) 
      { gmag += ds->grad[j]*ds->grad[j];
      }

      proj = 0;
      for (j = 0; j<ds->dim; j++) 
      { proj += ds->p[j]*ds->grad[j];
      }

      if (gmag>1e-30)
      { for (j = 0; j<ds->dim; j++)
        { ds->p[j] = ds->p[j] - 2*ds->grad[j]*proj/gmag;
        }      
      }
    }
  }

  it->proposals += 1;

  if (ds->pot_energy<=slice_point) 
  { for (j = 0; j<ds->dim; j++) ds->p[j] = - ds->p[j];
    /*fprintf(stderr," Accept %d\n");*/
  }
  else
  { mc_value_copy (ds->q, q_save, ds->dim);
    mc_value_copy (ds->grad, p_save, ds->dim);
    ds->pot_energy = old_pot;
    it->rejects += 1;
    /*fprintf(stderr," Reject\n");*/
  }
}
