/* MOL-UTIL.C - Molecular dynamics utility routines. */

/* Copyright (c) 1995-2003 by Radford M. Neal 
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
#include "mol.h"


/* COMPUTE SQUARED DISTANCE OF THE NEAREST IMAGES OF A PAIR OF MOLECULES. */

double squared_distance
( int D,		/* Number of dimensions */
  double len,		/* Length of each dimension */
  double *coords,	/* Coordinates of molecules */
  int i,		/* Indexes of molecules */
  int j
)
{
  double dx, dy, dz, d2, hlen;
  int iD, jD;

  hlen = len/2;
  iD = i*D;
  jD = j*D;

  dx = len * fmod (coords[iD]-coords[jD], 1.0);
  if (dx<0) dx = len + dx;
  if (dx>hlen) dx = len - dx; 
  d2 = dx*dx;
  if (D>1) 
  { dy = len * fmod (coords[iD+1]-coords[jD+1], 1.0); 
    if (dy<0) dy = len + dy;
    if (dy>hlen) dy = len - dy;
    d2 += dy*dy; 
    if (D>2) 
    { dz = len * fmod (coords[iD+2]-coords[jD+2], 1.0); 
      if (dz<0) dz = len + dz;
      if (dz>hlen) dz = len - dz; 
      d2 += dz*dz; 
    }
  }

  return d2;
}


/* COMPUTE ENERGY AND ITS GRADIENT FOR ONE PAIR OF PARTICLES.  The contribution
   from this pair to the energy is returned as the value of this procedure,
   and the gradients with respect to the coordinates (in [0,1) form) are stored
   in gradi and gradj, if they are non-null.  The gradient with respect to the 
   log of the length of a dimension is stored in gradl, if it is non-null. */

double energy_and_gradient
( double *coords,	/* Coordinates of molecules, in [0,1) scaling; 
			   log of length of a dimension at end, for NPT */
  double inv_temp,	/* Inverse temperature */
  mol_spec *ms,		/* Specifications of molecular system */
  int i,		/* Indexes of particles */
  int j,
  double *gradi,	/* Place to store gradient w.r.t. i's coordinates */
  double *gradj,	/* Place to store gradient w.r.t. j's coordinates */
  double *gradl		/* Place to store gradient w.r.t. log dimension length*/
)
{
  double w2, dx, dy, dz, d2, t2, t6, m, len, scale, max, hlen, dE;
  int D, nx, ny, nz, iD, jD, k;
  double energy;

  if (i==j) abort();

  D = ms->D;

  iD = i*D;
  jD = j*D;

  scale = ms->scale;

  if (inv_temp==0 || scale==0) 
  { for (k = 0; k<D; k++)
    { if (gradi) gradi[k] = 0;
      if (gradj) gradj[k] = 0;
    }
    return 0;
  }

  w2 = ms->width * ms->width;
  len = dlen(ms,coords);
  max = ms->max_pair_energy;
  hlen = len/2;

  /* Compute d2, the squared distance of the nearest images of the pair. */

  dx = len * fmod (coords[iD]-coords[jD], 1.0); 
  if (dx<0) dx = len + dx;
  nx = 0;
  if (dx>hlen) { dx = len - dx; nx = 1; }
  d2 = dx*dx;
  if (D>1) 
  { dy = len * fmod (coords[iD+1]-coords[jD+1], 1.0); 
    if (dy<0) dy = len + dy;
    ny = 0;
    if (dy>hlen) { dy = len - dy; ny = 1; }
    d2 += dy*dy; 
    if (D>2) 
    { dz = len * fmod (coords[iD+2]-coords[jD+2], 1.0); 
      if (dz<0) dz = len + dz;
      nz = 0; 
      if (dz>hlen) { dz = len - dz; nz = 1; } 
      d2 += dz*dz; 
    }
  }

  /* Compute the contribution to the energy and gradient of this pair. */

  t2 = w2 / d2;
  t6 = t2*t2*t2;
  dE = 4 * scale * t6*(t6-1);

  if (dE<=0 || dE<max-4*scale*d2/w2)
  {
    energy = inv_temp * dE;
    if (gradi || gradj)  
    { m = inv_temp * len * 48 * scale * t6*(t6-0.5) / d2;
    }
  }
  else
  { 
    energy = inv_temp * (max - 4*scale*d2/w2);
    if (gradi || gradj) 
    {  m = inv_temp * len * 8 * scale / w2;
    }
  }

  if (gradi)
  {
    gradi[0] = nx ? dx * m : -dx * m;

    if (D>1)
    { 
      gradi[1] = ny ? dy * m : -dy * m;

      if (D>2)
      { gradi[2] =  nz ? dz * m : -dz * m;
      }
    }
  }

  if (gradj)
  {
    gradj[0] = nx ? -dx * m : dx * m;

    if (D>1)
    { 
      gradj[1] = ny ? -dy * m : dy * m;

      if (D>2)
      { gradj[2] =  nz ? -dz * m : dz * m;
      }
    }
  }

  if (ms->len_pres<0 && gradl)
  { *gradl = fabs(coords[ms->N*ms->D])>100 ? 0 :- m * d2 / len;
  }

  return energy;
}
