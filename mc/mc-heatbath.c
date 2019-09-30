/* MC-HEATBATH.C - Procedures for performing heatbath updates. */

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


/* PERFORM HEATBATH UPDATE OF MOMENTUM. */

void mc_heatbath
( mc_dynamic_state *ds,	/* State to update */
  float temperature,	/* Temperature of heat bath */
  float decay		/* Decay factor for existing momentum */
)
{ 
  double std_dev;
  int j;

  std_dev = sqrt (temperature * (1-decay*decay));

  if (decay==0)
  { for (j = 0; j<ds->dim; j++)
    { ds->p[j] = std_dev * rand_gaussian();
    }
  }
  else
  { for (j = 0; j<ds->dim; j++)
    { ds->p[j] = decay * ds->p[j] + std_dev * rand_gaussian();
    }
  }

  ds->know_kinetic = 0;
}


/* PERFORM RADIAL HEATBATH UPDATE OF MOMENTUM. */

void mc_radial_heatbath
( mc_dynamic_state *ds,	/* State to update */
  float temperature	/* Temperature of heat bath */
)
{ 
  double new_kinetic, F;
  int j;

  if (!ds->know_kinetic)
  { ds->kinetic_energy = mc_kinetic_energy(ds);
    ds->know_kinetic = 1;
  }

  new_kinetic = temperature * rand_gamma(ds->dim/2);

  F = sqrt(new_kinetic/ds->kinetic_energy);

  for (j = 0; j<ds->dim; j++)
  { ds->p[j] *= F;
  }

  ds->kinetic_energy = new_kinetic;
}


/* MIX MOMENTUM WITHOUT CHANGING TOTAL. */

void mc_mix_momentum
( mc_dynamic_state *ds,	/* State to update */
  float mix		/* Amount to mix by */
)
{
  double r, f;
  int i;

  if (!ds->know_kinetic)
  { ds->kinetic_energy = mc_kinetic_energy(ds);
    ds->know_kinetic = 1;
  }

  r = sqrt(2*ds->kinetic_energy);

  for (i = 0; i<ds->dim; i++)
  { ds->p[i] += mix * (r/sqrt(ds->dim)) * rand_gaussian();
  } 

  f = r / sqrt(2*mc_kinetic_energy(ds));

  for (i = 0; i<ds->dim; i++)
  { ds->p[i] *= f;
  } 
}
