/* MC-UTIL.C - Utility procedures for Monte Carlo simulation. */

/* Copyright (c) 1995-2019 by Radford M. Neal 
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


/* POINTER TO APPLICATION'S PROCEDURE FOR SETTING CUSTOM RANGES.  Zero if
   not supplied by application.  The range is specified by the second 
   argument, which is -2 or below, for 'A' to 'Z'.  The value returned is
   zero if the range specification is unknown to the application.  */

int (*mc_app_set_range_ptr) (mc_dynamic_state *, int *, int *);


/* SET REQUIRED RECORD SIZES. */

void mc_record_sizes
( log_gobbled *logg
)
{
  mc_app_record_sizes(logg);

  logg->req_size['r'] = sizeof (rand_state);
  logg->req_size['i'] = sizeof (mc_iter);
  logg->req_size['o'] = sizeof (mc_ops);
  logg->req_size['t'] = sizeof (mc_traj);
  logg->req_size['b'] = sizeof (mc_temp_state);
  logg->req_size['m'] = sizeof (mc_temp_sched);
  logg->req_size['h'] = sizeof (mc_therm_state);
}


/* COMPUTE KINETIC ENERGY OF MOMENTUM VARIABLES.  Returns zero if no momentum
   variables exist. */

double mc_kinetic_energy
( mc_dynamic_state *ds
)
{
  double s;
  int j;

  if (ds->p==0) return 0;

  s = 0;

  for (j = 0; j<ds->dim; j++)
  { s += ds->p[j] * ds->p[j];
  }

  return s / 2;
}


/* ENSURE THERMOSTAT IS PRESENT.  Currently disabled. */

#if 0

void mc_therm_present
( mc_dynamic_state *ds	/* State to update */
)
{
  static mc_therm_state thm;

  if (ds->therm_state==0)
  {
    ds->therm_state = &thm;

    ds->therm_state->tq = 1;
    ds->therm_state->tp = 0;
  }
}

#endif


/* ENSURE TEMPERING RECORD IS PRESENT. */

void mc_temp_present
( mc_dynamic_state *ds,	/* State to update */
  mc_temp_sched *sch	/* Tempering schedule */
)
{
  static mc_temp_state ts;

  if (ds->temp_state==0)
  {
    ds->temp_state = &ts;

    ds->temp_state->inv_temp = 1;
    ds->temp_state->temp_dir = 1;

    ds->temp_index = mc_temp_index (sch, ds->temp_state->inv_temp);
  }
}


/* FIND INDEX FOR SIMULATED TEMPERING INVERSE TEMPERATURE. */

int mc_temp_index
( mc_temp_sched *sch,		/* Schedule of inverse temperatures */
  float inv_temp		/* Inverse temperature to find */
)
{ 
  int i;

  if (sch==0)
  { fprintf(stderr,"No tempering schedule has been specified\n");
    exit(1);
  }

  for (i = 0; sch->sched[i].inv_temp!=inv_temp; i++)
  { if (i==Max_temps-1) abort();
  }

  return i;
}


/* FIND [0,1] VARIABLE FOR ACCEPT/REJECT OR SLICE LEVEL.  Updates or 
   records in slevel as required. */

double mc_slevel
( mc_dynamic_state *ds		/* Dynamical state structure */
)
{ double a;
  if (ds->slevel.random==-1 || ds->slevel.random==2)
  { a = rand_uniform();
  }
  else
  { a = ds->slevel.value;
    if (ds->slevel.random!=0)
    { a = a + ds->slevel.random * (2*rand_uniform()-1);
    }
    a = a + ds->slevel.move;
    while (a > 1) a -= 2;
    while (a < -1) a += 2;
  }
  ds->slevel.value = a;
  return a>=0 ? a : -a;
}


/* FIND SLICE LEVEL INCREMENT. */

double mc_slice_inc
( mc_dynamic_state *ds		/* Dynamical state structure */
)
{
  double s = mc_slevel(ds);
  if (s < 1e-30) s = 1e-30;
  return -log(s);
}


/* FIND DIFFERENCE IN ENERGY FOR TEMPERING CHANGE.  Computes the difference
   in energy resulting from a change in the tempering index by one step
   up or down. */

double mc_energy_diff
( mc_dynamic_state *ds,		/* Dynamical state structure */
  mc_temp_sched *sch,		/* Tempering schedule */
  int dir			/* Direction of change */
)
{ 
  double ed, e2;

  if (ds->temp_state==0 
   || ds->temp_state->inv_temp==1 && dir>0
   || ds->temp_index==0 && dir<0)
  { abort();
  }

  if (!ds->know_pot) 
  { mc_app_energy(ds,1,1,&ds->pot_energy,0);
    ds->know_pot = 1;
  }
  ds->temp_index += dir;
  ds->temp_state->inv_temp = sch->sched[ds->temp_index].inv_temp;
  mc_app_energy(ds,1,1,&e2,0);
  ds->temp_index -= dir;
  ds->temp_state->inv_temp = sch->sched[ds->temp_index].inv_temp;
  ed = e2 - ds->pot_energy;

  return ed;
}


/* COPY STATE VARIABLES. */

void mc_value_copy
( mc_value *dest,	/* Place to copy to */
  mc_value *src,	/* Place to copy from */
  int n			/* Number of values to copy */
)
{
  while (n>0)
  { *dest++ = *src++;
    n -= 1;
  }
}


/* SET RANGE OF COMPONENTS FOR UPDATES. */

void mc_set_range
( mc_dynamic_state *ds,
  int *first,
  int *last,
  int r_update
)
{
  if (*first>=0)
  { if (*last>ds->dim-1)
    { fprintf (stderr, "Last (%d) is beyond end (%d) in coordinate range\n",
               *last, ds->dim-1);
      exit(1);
    }
  }
  else if (*first==-1)
  { *first = 0;
    *last = ds->dim-1;
  }
  else 
  { if (mc_app_set_range_ptr==0 || !(*mc_app_set_range_ptr) (ds, first, last))
    { fprintf (stderr, "Unknown application-specified coordinate range (%c)\n",
                       'A' + (-2-*first));
      exit(1);
    }
  }

  if (*first>*last+1)
  { fprintf (stderr, 
     "First (%d) greater than last (%d) + 1 in coordinate range\n",
      *first, *last);
    exit(1);
  }

  if (r_update)
  { *first = *last = *first + (int)(rand_uniform()*(*last-*first+1));
  }
}
