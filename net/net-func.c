/* NET-FUNC.C - Routine for calculating the function defined by a network. */

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
#include "net.h"


/* This module calculates the values of the output units in a network, given 
   values for the input units.  The values of hidden units are calculated
   along the way.  There are facilities for starting the calculation on the 
   assumption the values are already known up to some layer, as would be 
   the case if the weights into earlier layers have not changed since the
   last calculation. 
*/


static void bias_values (net_value *, int, net_param *);

static void add_connections (net_value *, int, net_value *, int, 
                             net_param *, net_param *);


/* EVALUATE NETWORK FUNCTION FOR GIVEN INPUTS.  The inputs are taken from
   the net_values structure passed.  When 'start' is greater than zero, the
   correct unit values for that number of hidden layers are assumed to be
   already present in the net_values structure. */

void net_func 
( net_values *v,	/* Place to get inputs and store outputs */
  int start,		/* Number of hidden layers with known values */
  net_arch *a,		/* Network architecture */
  net_params *w		/* Network parameters */
)
{
  net_value *vh, *sh;
  int l, j;

  /* Compute values for successive hidden layers. */

  for (l = start; l<a->N_layers; l++)
  {
    sh = v->s[l];
    vh = v->h[l];

    bias_values (sh, a->N_hidden[l], a->has_bh[l] ? w->bh[l] : 0);

    if (a->has_ih[l])
    { add_connections (sh, a->N_hidden[l], v->i, a->N_inputs, 
                       w->ih[l], a->has_ti ? w->ti : 0);
    }

    if (l>0 && a->has_hh[l-1])
    { add_connections (sh, a->N_hidden[l], v->h[l-1], a->N_hidden[l-1],
                       w->hh[l-1], a->has_th[l-1] ? w->th[l-1] : 0);
    }

    /* Put values through 'tanh'. */

    for (j = 0; j<a->N_hidden[l]; j++)
    { vh[j] = tanh(sh[j]);
    }
  }

  /* Compute values for the outputs. */

  bias_values (v->o, a->N_outputs, a->has_bo ? w->bo : 0);

  if (a->has_io)
  { add_connections (v->o, a->N_outputs, v->i, a->N_inputs,
                     w->io, a->has_ti ? w->ti : 0);
  }

  for (l = 0; l<a->N_layers; l++)
  {
    if (a->has_ho[l])
    { add_connections (v->o, a->N_outputs, v->h[l], a->N_hidden[l], 
                       w->ho[l], a->has_th[l] ? w->th[l] : 0);
    }
  }
}


/* SET UNIT VALUES TO BIASES.  Just zeros them if there are no biases. */

static void bias_values
( net_value *v,		/* Array of unit values to set */
  int n,		/* Number of units */
  net_param *b		/* Biases, null if none */
)
{
  int j;

  if (b!=0)
  {
    for (j = 0; j<n; j++) v[j] = *b++;
  }
  else
  {
    for (j = 0; j<n; j++) v[j] = 0;
  }
}


/* ADD CONTRIBUTION FROM ONE GROUP OF CONNECTIONS.  Adds the weighted input
   due to connections from one source layer to the current unit values for
   the destination layer. */

static void add_connections
( net_value *s,		/* Summed input for destination units to add to */
  int nd,		/* Number of destination units */
  net_value *v,		/* Values for source units */
  int ns,		/* Number of source units */
  net_param *w,		/* Connection weights */
  net_param *t		/* Offsets to add to source unit values */
)
{
  net_value tv;
  int i, j;

  if (t!=0)
  {
    for (i = 0; i<ns; i++)
    { tv = v[i] + *t++;
      j = 0;
      do { s[j] += *w++ * tv; j += 1; } while (j<nd);
    }
  }
  else
  {
    for (i = 0; i<ns; i++)
    { tv = v[i];
      j = 0;
      do { s[j] += *w++ * tv; j += 1; } while (j<nd);
    }
  }
}
