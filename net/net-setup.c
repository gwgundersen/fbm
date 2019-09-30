/* NET-SETUP.C - Procedures for setting up network data structures. */

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


/* This module contains routines for setting up the structures containing
   pointers to arrays that are used to describe the parameters, hyperparameters,
   and unit values for networks.  

   The scheme used allows all quantites of a given type to be stored in
   a contiguous block of storage.  The pointer structure allows the various 
   sections of this block to be easily accessed. 

   Sub-groups of parameters or hyperparameters may also be located via
   a group index, using the net_setup_hyper_group and net_setup_param_group
   procedures.
*/


/* RETURN NUMBER OF HYPERPARAMETERS.  Returns the number of 'sigma' values
   for a network with the given architecture. */
   
int net_setup_sigma_count
( net_arch *a,		/* Network architecture */
  net_flags *flgs,	/* Network flags, null if none */
  model_specification *m /* Data model */
)
{ 
  int count;
  int l;

  count = 0;

  if (a->has_ti) count += 1;
  
  for (l = 0; l<a->N_layers; l++)
  { 
    if (l>0 && a->has_hh[l-1]) count += 1 + a->N_hidden[l-1];

    if (a->has_ih[l]) 
    { count += 1 + not_omitted(flgs?flgs->omit:0,a->N_inputs,1<<(l+1));
    }
    if (a->has_bh[l]) count += 1;
    if (a->has_ah[l]) count += a->N_hidden[l];
    if (a->has_th[l]) count += 1;
    if (a->has_ho[l]) count += 1 + a->N_hidden[l];
  }

  if (a->has_io) 
  { count += 1 + not_omitted(flgs?flgs->omit:0,a->N_inputs,1);
  }
  if (a->has_bo) count += 1;
  if (a->has_ao) count += a->N_outputs;

  if (m!=0 && m->type=='R') count += 1 + a->N_outputs;
  
  return count;
}


/* RETURN NUMBER OF NETWORK PARAMETERS.  Returns the number of weights,
   biases, and offsets in a network with the given architecture. */

int net_setup_param_count
( net_arch *a,		/* Network architecture */
  net_flags *flgs	/* Network flags, null if none */
)
{
  int count;
  int l;

  count = 0;
 
  if (a->has_ti) count += a->N_inputs;

  for (l = 0; l<a->N_layers; l++)
  {
    if (l>0 && a->has_hh[l-1]) count += a->N_hidden[l-1]*a->N_hidden[l];

    if (a->has_ih[l]) 
    { count += not_omitted(flgs?flgs->omit:0,a->N_inputs,1<<(l+1))
                 * a->N_hidden[l];
    }

    if (a->has_bh[l]) count += a->N_hidden[l];
    if (a->has_th[l]) count += a->N_hidden[l];
    if (a->has_ho[l]) count += a->N_hidden[l]*a->N_outputs;
  }

  if (a->has_io) 
  { count += not_omitted(flgs?flgs->omit:0,a->N_inputs,1)
              * a->N_outputs;
  }

  if (a->has_bo) count += a->N_outputs;
  
  return count;
}


/* RETURN NUMBER OF UNIT-RELATED VALUES IN NETWORK.  Returns the number 
   of unit-related values for a network with the given architecture. */

int net_setup_value_count 
( net_arch *a		/* Network architecture */
)
{ 
  int count;
  int l;

  count = a->N_inputs + a->N_outputs;

  for (l = 0; l<a->N_layers; l++)
  { count += 2 * a->N_hidden[l];
  }

  return count;
}


/* SET UP POINTERS TO HYPERPARAMETERS.  Sets the pointers in the net_sigmas
   structure to point the appropriate places in a block of sigma values.
   Pointers associated with parts of the network that don't exist are set to 
   zero.  

   The size of the storage block, as returned by the net_setup_sigma_count 
   procedure, must be stored by the caller in the total_sigmas field, and
   the caller must also set the sigma_block field to point to a block of this
   many net_sigma values. */

void net_setup_sigma_pointers
( net_sigmas *s,	/* Structure to set up pointers in */
  net_arch *a,		/* Network architecture */
  net_flags *flgs,	/* Network flags, null if none */
  model_specification *m /* Data model */
)
{ 
  net_sigma *b;
  int l;

  b = s->sigma_block;

  s->ti_cm = a->has_ti ? b++ : 0;

  for (l = 0; l<a->N_layers; l++)
  {
    if (l>0)
    { s->hh_cm[l-1] = a->has_hh[l-1] ? b++ : 0;
      s->hh[l-1] = 0;
      if (a->has_hh[l-1]) 
      { s->hh[l-1] = b;
        b += a->N_hidden[l-1];
      }
    }

    s->ih_cm[l] = a->has_ih[l] ? b++ : 0;
    s->ih[l] = 0;
    if (a->has_ih[l]) 
    { s->ih[l] = b;
      b += not_omitted(flgs?flgs->omit:0,a->N_inputs,1<<(l+1));
    }

    s->bh_cm[l] = a->has_bh[l] ? b++ : 0;

    s->ah[l] = 0;
    if (a->has_ah[l])
    { s->ah[l] = b;
      b += a->N_hidden[l];
    }
  
    s->th_cm[l] = a->has_th[l] ? b++ : 0;
  }

  for (l = a->N_layers-1; l>=0; l--)
  { s->ho_cm[l] = a->has_ho[l] ? b++ : 0;
    s->ho[l] = 0;
    if (a->has_ho[l]) 
    { s->ho[l] = b;
      b += a->N_hidden[l];
    }
  }

  s->io_cm = a->has_io ? b++ : 0;
  s->io = 0;
  if (a->has_io) 
  { s->io = b;
    b += not_omitted(flgs?flgs->omit:0,a->N_inputs,1);
  }

  s->bo_cm = a->has_bo ? b++ : 0;

  s->ao = 0;
  if (a->has_ao)
  { s->ao = b;
    b += a->N_outputs;
  }

  s->noise_cm = m!=0 && m->type=='R' ? b++ : 0;
  s->noise = 0;
  if (m!=0 && m->type=='R') 
  { s->noise = b;
    b += a->N_outputs;
  }
}


/* SET UP POINTERS TO NETWORK PARAMETERS.  Sets the pointers in the net_params
   structure to point the appropriate places in a block of parameter values.
   Pointers associated with parts of the network that don't exist are set to 
   zero.  

   The size of the storage block, as returned by the net_setup_param_count
   procedure, must be stored by the caller in the total_params field, and
   the caller must also set the param_block field to point to a block of this 
   many net_param values. */

void net_setup_param_pointers
( net_params *w,	/* Structure to set up pointers in */
  net_arch *a,		/* Network architecture */
  net_flags *flgs	/* Network flags, null if none */
)
{
  net_param *b;
  int l;

  b = w->param_block;

  w->ti = 0;
  if (a->has_ti)
  { w->ti = b;
    b += a->N_inputs;
  }

  for (l = 0; l<a->N_layers; l++)
  {
    if (l>0)
    { w->hh[l-1] = 0;
      if (a->has_hh[l-1]) 
      { w->hh[l-1] = b;
        b += a->N_hidden[l-1]*a->N_hidden[l];
      }
    }

    w->ih[l] = 0;
    if (a->has_ih[l]) 
    { w->ih[l] = b;
      b += not_omitted(flgs?flgs->omit:0,a->N_inputs,1<<(l+1))
            * a->N_hidden[l];
    }
  
    w->bh[l] = 0;
    if (a->has_bh[l])
    { w->bh[l] = b;
      b += a->N_hidden[l];
    }
  
    w->th[l] = 0;
    if (a->has_th[l])
    { w->th[l] = b;
      b += a->N_hidden[l];
    }
  }

  for (l = a->N_layers-1; l>=0; l--)
  { w->ho[l] = 0;
    if (a->has_ho[l]) 
    { w->ho[l] = b;
      b += a->N_hidden[l]*a->N_outputs;
    }
  }

  w->io = 0;
  if (a->has_io) 
  { w->io = b;
    b += not_omitted(flgs?flgs->omit:0,a->N_inputs,1) * a->N_outputs;
  }

  w->bo = 0;
  if (a->has_bo)
  { w->bo = b;
    b += a->N_outputs;
  }
}


/* SET UP POINTERS TO UNIT VALUES.  Sets the pointers in the net_values
   structure to point the appropriate places in the block of net_value
   values passed.  The size of this block must be as indicated by the
   net_setup_value_count procedure, which count should be stored by
   the caller in the total_values field.. */

void net_setup_value_pointers
( net_values *v,	/* Structure to set up pointers in */
  net_value *b,		/* Block of 'value' values */
  net_arch *a		/* Network architecture */
)
{
  int l;

  v->i = b;
  b += a->N_inputs;

  for (l = 0; l<a->N_layers; l++)
  { v->h[l] = b;
    b += a->N_hidden[l];
    v->s[l] = b;
    b += a->N_hidden[l];
  }

  v->o = b;
  b += a->N_outputs;
}


/* LOCATE GROUP OF HYPERPARAMETERS.  Finds the offset and number of 
   hyperparameters in a group, plus whether it's a group of adjustments.  
   Groups are identified by integers from one on up.  Zero is returned 
   if the group index is too big. 

   Noise sigmas are not considered hyperparameters for this procedure. */

int net_setup_hyper_group
( net_arch *a,		/* Network architecture */
  net_flags *flgs,	/* Network flags, null if none */
  int grp,		/* Index of group */
  int *offset,		/* Set to offset of group within block */
  int *number,		/* Set to number of items in group */
  int *adj		/* Set to whether this is a group of adjustments */
)
{ 
  int i, l;

  *adj = 0;

  if (grp<1) return 0;

  i = 0;

  if (a->has_ti) 
  { *offset = i; 
    i += 1; 
    if (--grp==0) goto done;
  }

  for (l = 0; l<a->N_layers; l++)
  {
    if (l>0)
    { if (a->has_hh[l-1]) 
      { *offset = i; 
        i += 1 + a->N_hidden[l-1]; 
        if (--grp==0) goto done;
      }
    }

    if (a->has_ih[l]) 
    { *offset = i; 
      i += 1 + not_omitted(flgs?flgs->omit:0,a->N_inputs,1<<(l+1)); 
      if (--grp==0) goto done;
    }

    if (a->has_bh[l]) 
    { *offset = i; 
      i += 1; 
      if (--grp==0) goto done;
    }

    if (a->has_ah[l]) 
    { *offset = i-1;  /* Pretend that there is a high-level one before */ 
      i += a->N_hidden[l]; 
      if (--grp==0) { *adj = 1; goto done; }
    }

    if (a->has_th[l]) 
    { *offset = i; 
      i += 1; 
      if (--grp==0) goto done;
    }
  }

  for (l = a->N_layers-1; l>=0; l--)
  { if (a->has_ho[l]) 
    { *offset = i; 
      i += 1 + a->N_hidden[l]; 
      if (--grp==0) goto done;   
    }
  }

  if (a->has_io) 
  { *offset = i;   
    i += 1 + not_omitted(flgs?flgs->omit:0,a->N_inputs,1); 
    if (--grp==0) goto done;
  }

  if (a->has_bo) 
  { *offset = i; 
    i += 1; 
    if (--grp==0) goto done;  
  }

  if (a->has_ao) 
  { *offset = i-1;  /* Pretend that there is a high-level one before */
    i += a->N_outputs; 
    if (--grp==0) { *adj = 1; goto done; }
  }

  return 0;

done:
  *number = i-*offset;
  return 1;
}


/* LOCATE GROUP OF PARAMETERS.  Finds the offset, number, and dimension of 
   for a group of parameters.  Groups are identified by integers 
   from one on up.  Zero is returned if the group index is too big,
   or refers to an adjustment group (having no parameters in it). */

int net_setup_param_group
( net_arch *a,		/* Network architecture */
  net_flags *flgs,	/* Network flags, null if none */
  int grp,		/* Index of group */
  int *offset,		/* Set to offset of group within block */
  int *number,		/* Set to number of items in group */
  int *source		/* Set to number of source units associated with group,
                             or to zero if it's a one-dimensional group */
)
{ 
  int i, l;

  if (grp<1) return 0;

  i = 0;

  if (a->has_ti) 
  { *offset = i; 
    *source = 0;
    i += a->N_inputs; 
    if (--grp==0) goto done;
  }

  for (l = 0; l<a->N_layers; l++)
  {
    if (l>0)
    { if (a->has_hh[l-1]) 
      { *offset = i; 
        *source = a->N_hidden[l-1];
        i += a->N_hidden[l-1]*a->N_hidden[l]; 
        if (--grp==0) goto done;
      }
    }

    if (a->has_ih[l]) 
    { *offset = i; 
      *source = not_omitted(flgs?flgs->omit:0,a->N_inputs,1<<(l+1));
      i += *source * a->N_hidden[l]; 
      if (--grp==0) goto done;
    }

    if (a->has_bh[l]) 
    { *offset = i; 
      *source = 0;
      i += a->N_hidden[l]; 
      if (--grp==0) goto done;
    }

    if (a->has_ah[l])
    { if (--grp==0) return 0;
    }

    if (a->has_th[l]) 
    { *offset = i; 
      *source = 0;
      i += a->N_hidden[l]; 
      if (--grp==0) goto done;
    }
  }

  for (l = a->N_layers-1; l>=0; l--)
  { if (a->has_ho[l]) 
    { *offset = i; 
      *source = a->N_hidden[l];
      i += a->N_hidden[l]*a->N_outputs; 
      if (--grp==0) goto done;
    }
  }

  if (a->has_io) 
  { *offset = i;    
    *source = not_omitted(flgs?flgs->omit:0,a->N_inputs,1);
    i += *source * a->N_outputs; 
    if (--grp==0) goto done;
  }

  if (a->has_bo) 
  { *offset = i; 
    *source = 0;
    i += a->N_outputs; 
    if (--grp==0) goto done;
  }

  if (a->has_ao)
  { if (--grp==0) return 0;
  }

  return 0;

done:
  *number = i-*offset;
  return 1;
}
