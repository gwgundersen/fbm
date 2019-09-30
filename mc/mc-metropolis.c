/* MC-METROPOLIS.C - Procedures for performing Metropolis-style updates. */

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


/* PERFORM METROPOLIS UPDATE ON ALL COMPONENTS AT ONCE. */

void mc_metropolis
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  mc_value *q_save,	/* Place to save old q values */
  int b_accept		/* Use Barker/Boltzmann acceptance probability? */
)
{
  double old_energy, sf, U, a;
  int k;

  if (!ds->know_pot)
  { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
    ds->know_pot = 1;
  }

  old_energy = ds->pot_energy;

  sf = it->stepsize_factor;

  mc_value_copy (q_save, ds->q, ds->dim);

  for (k = 0; k<ds->dim; k++) 
  { ds->q[k] += sf * ds->stepsize[k] * rand_gaussian();
  }
  
  mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);

  it->proposals += 1;
  it->delta = ds->pot_energy - old_energy;

  U = rand_uniform(); /* Do every time to keep in sync for coupling purposes */

  a = exp(-it->delta/it->temperature);
  if (b_accept) 
  { a = 1/(1+1/a);
  }

  if (U<a)
  { 
    it->move_point = 1;
    ds->know_grad = 0;
  }
  else
  { 
    it->rejects += 1;
    it->move_point = 0;

    ds->pot_energy = old_energy;

    mc_value_copy (ds->q, q_save, ds->dim);
  }
}


/* PERFORM RANDOM-GRID METROPOLIS UPDATE ON ALL COMPONENTS. */

void mc_rgrid_met
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  mc_value *q_save,	/* Place to save old q values */
  int b_accept		/* Use Barker/Boltzmann acceptance probability? */
)
{
  double old_energy, sf, U, a, ss;
  int k;

  if (!ds->know_pot)
  { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
    ds->know_pot = 1;
  }

  old_energy = ds->pot_energy;

  sf = it->stepsize_factor;

  mc_value_copy (q_save, ds->q, ds->dim);

  for (k = 0; k<ds->dim; k++) 
  { U = rand_uniopen() - 0.5;
    ss = sf * ds->stepsize[k];
    ds->q[k] = (2*ss) * (U + floor (0.5 + ds->q[k]/(2*ss) - U));
  }
  
  mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);

  it->proposals += 1;
  it->delta = ds->pot_energy - old_energy;

  U = rand_uniform(); /* Do every time to keep in sync for coupling purposes */

  a = exp(-it->delta/it->temperature);
  if (b_accept) 
  { a = 1/(1+1/a);
  }

  if (U<a)
  { 
    it->move_point = 1;
    ds->know_grad = 0;
  }
  else
  { 
    it->rejects += 1;
    it->move_point = 0;

    ds->pot_energy = old_energy;

    mc_value_copy (ds->q, q_save, ds->dim);
  }
}


/* PERFORM METROPOLIS UPDATE ON ONE COMPONENT AT A TIME. */

void mc_met_1
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int firsti,		/* Index of first component to update (-1 for all) */
  int lasti,		/* Index of last component to update */
  int b_accept,		/* Use Barker/Boltzmann acceptance probability? */
  int r_update		/* Update just one component at random? */
)
{
  double old_energy, qsave, sf, U, a;
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

  sf = it->stepsize_factor;

  for (k = firsti; k<=lasti; k++)
  {
    if (!ds->know_pot)
    { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
      ds->know_pot = 1;
    }
  
    old_energy = ds->pot_energy;
  
    qsave = ds->q[k];
  
    ds->q[k] += sf * ds->stepsize[k] * rand_gaussian();
    
    mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
  
    it->proposals += 1;
    it->delta = ds->pot_energy - old_energy;
  
    U = rand_uniform(); /* Do every time to keep in sync for coupling purposes*/

    a = exp(-it->delta/it->temperature);
    if (b_accept) 
    { a = 1/(1+1/a);
    }

    if (U<a)
    { 
      it->move_point = 1;
      ds->know_grad = 0;
    }
    else
    { 
      it->rejects += 1;
      it->move_point = 0;
  
      ds->pot_energy = old_energy;
  
      ds->q[k] = qsave;
    }
  }
}


/* PERFORM RANDOM-GRID METROPOLIS UPDATE ON ONE COMPONENT AT A TIME. */

void mc_rgrid_met_1
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int firsti,		/* Index of first component to update (-1 for all) */
  int lasti,		/* Index of last component to update */
  int b_accept,		/* Use Barker/Boltzmann acceptance probability? */
  int r_update		/* Update just one component at random? */
)
{
  double old_energy, qsave, sf, ss, U, a;
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

  sf = it->stepsize_factor;

  for (k = firsti; k<=lasti; k++)
  {
    if (!ds->know_pot)
    { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
      ds->know_pot = 1;
    }
  
    old_energy = ds->pot_energy;
  
    qsave = ds->q[k];
  
    U = rand_uniopen() - 0.5;
    ss = sf * ds->stepsize[k];
    ds->q[k] = (2*ss) * (U + floor (0.5 + ds->q[k]/(2*ss) - U));
    
    mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
  
    it->proposals += 1;
    it->delta = ds->pot_energy - old_energy;
  
    U = rand_uniform(); /* Do every time to keep in sync for coupling purposes*/

    a = exp(-it->delta/it->temperature);
    if (b_accept) 
    { a = 1/(1+1/a);
    }

    if (U<a)
    { 
      it->move_point = 1;
      ds->know_grad = 0;
    }
    else
    { 
      it->rejects += 1;
      it->move_point = 0;
  
      ds->pot_energy = old_energy;
  
      ds->q[k] = qsave;
    }
  }
}


/* PERFORM GAUSSIAN GIBBS SAMPLING UPDATES.  The components updated are ASSUMED,
   without checking, to have Gaussian conditional distributions, which are 
   found by evaluating the energy at points -1, 0, +1. */

void mc_gaussian_gibbs
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int firsti,		/* Index of first component to update (-1 for all) */
  int lasti,		/* Index of last component to update */
  int r_update		/* Update just one component at random? */
)
{
  double Eminus, Ezero, Eplus;
  double mean, tau;
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
    ds->q[k] = -1;
    mc_app_energy (ds, 1, 1, &Eminus, 0);

    ds->q[k] = 0;
    mc_app_energy (ds, 1, 1, &Ezero, 0);

    ds->q[k] = 1;
    mc_app_energy (ds, 1, 1, &Eplus, 0);

    tau = Eminus - 2*Ezero + Eplus;
    if (tau<=0)
    { fprintf(stderr,"Variance in gaussian_gibbs is not positive!\n");
      exit(1);
    }

    mean = (Eminus - Eplus) / (2*tau);

    ds->q[k] = mean + rand_gaussian() / sqrt(tau);
  }

  ds->know_pot = 0;
  ds->know_grad = 0;
}
