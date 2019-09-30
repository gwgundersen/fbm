/* NET-BACK.C - Routine for backpropagating the "error" through the network. */

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
#include "net.h"


/* This module finds the derivative of the "error" for a particular case
   with respect to the values of the hidden and input units, and the with
   respect to the summed input for the hidden units, given the derivative 
   with respect to the output units.  There are facilities for calculating 
   these derivatives only back to a certain layer.
*/

#define sqrt_2 1.4142135623730950488

static void zero_derivatives (net_value *, int),
            sum_derivatives  (net_value *, int, net_value *, int, net_param *, 
                              char *, int);


/* BACKPROPAGATE ERROR DERIVATIVES.  The first argument must contain the 
   values of all the units in the network.  The second must contain the
   derivatives of the "error" with respect to the output units (in g->o).
   These derivatives are backpropagated to give the derivatives of the
   error with respect to the other units in the network, and with respect
   to the summed input into the hidden units..  This is done back to hidden 
   layer 'start', or back all the way to the inputs if 'start' is -1. */

void net_back
( net_values *v,	/* Values for units in network */
  net_values *d,	/* Place to get output derivatives, and store others */
  int start,		/* Earliest layer to find derivatives for */
  net_arch *a,		/* Network architecture */
  net_flags *flgs,	/* Network flags, null if none */
  net_params *w		/* Network parameters */
)
{
  int l, i;

  /* Backpropagate through hidden layers. */

  for (l = a->N_layers-1; l>=0 && l>=start; l--)
  { 
    zero_derivatives (d->h[l], a->N_hidden[l]);
    
    if (a->has_ho[l])
    { sum_derivatives (d->o, a->N_outputs, d->h[l], a->N_hidden[l], 
                       w->ho[l], (char *) 0, 0);
    }

    if (l<a->N_layers-1 && a->has_hh[l])
    { sum_derivatives (d->s[l+1], a->N_hidden[l+1], d->h[l], a->N_hidden[l], 
                       w->hh[l], (char *) 0, 0);
    }

    switch (flgs==0 ? Tanh_type : flgs->layer_type[l])
    { case Tanh_type:
      { for (i = 0; i<a->N_hidden[l]; i++)
        { d->s[l][i] = (1 - v->h[l][i]*v->h[l][i]) * d->h[l][i];
        }
        break;
      }
      case Sin_type:
      { for (i = 0; i<a->N_hidden[l]; i++)
        { d->s[l][i] = 2 * cos(v->s[l][i]*sqrt_2) * d->h[l][i];
        }
        break;
      }
      case Identity_type: 
      { for (i = 0; i<a->N_hidden[l]; i++)
        { d->s[l][i] = d->h[l][i];
        }
        break;
      }
      default: abort();
    }
  }

  /* Backpropagate to input layer. */

  if (start<0)
  {
    zero_derivatives (d->i, a->N_inputs);

    if (a->has_io)
    { sum_derivatives (d->o, a->N_outputs, d->i, a->N_inputs, w->io,
                       flgs ? flgs->omit : 0, 1);
    }
 
    for (l = 0; l<a->N_layers; l++)
    { if (a->has_ih[l])
      { sum_derivatives (d->s[l], a->N_hidden[l], d->i, a->N_inputs, w->ih[l],
                         flgs ? flgs->omit : 0, 1<<(l+1));
      }
    }
  }
}


/* ZERO DERIVATIVES.  Sets the derivatives with respect to a set of source
   units to zero. */

static void zero_derivatives
( net_value *d,		/* Derivatives with respect to units to zero */
  int n			/* Number of units */
)
{
  int i;
  
  for (i = 0; i<n; i++)
  { d[i] = 0;
  }
}


/* SUM UP CONTRIBUTIONS TO THE DERIVATIVES FROM ONE GROUP OF CONNECTIONS.  Adds 
   the weighted sum of derivatives due to connections from source units to 
   a given destination layer to the totals for the source layer. */

static void sum_derivatives
( net_value *dd,	/* Derivatives with respect to destination units */
  int nd,		/* Number of destination units */
  net_value *ds,	/* Derivatives w.r.t. source units to add to */
  int ns,		/* Number of source units */
  net_param *w,		/* Connection weights */
  char *omit,		/* Omit flags, null if not present */
  int b			/* Bit to look at in omit flags */
)
{
  net_value tv;
  int i, j, k;

  if (omit==0)
  {
    if (nd==1)
    { 
      for (i = 0; i<ns; i++)
      { ds[i] += *w++ * dd[0];
      }
    }
    else
    {
      for (i = 0; i<ns; i++)
      { tv = *w++ * dd[0];
        j = 1;
        do { tv += *w++ * dd[j]; j += 1; } while (j<nd);
        ds[i] += tv;
      }
    }
  }
  else
  {
    if (nd==1)
    { k = 0;
      for (i = 0; i<ns; i++)
      { if ((omit[i]&b)==0)
        { ds[k] += *w++ * dd[0];
          k += 1;
        }
      }
    }
    else
    {
      k = 0;
      for (i = 0; i<ns; i++)
      { if ((omit[i]&b)==0)
        { tv = *w++ * dd[0];
          j = 1;
          do { tv += *w++ * dd[j]; j += 1; } while (j<nd);
          ds[k] += tv;
          k += 1;
        }
      }
    }
  }

}
