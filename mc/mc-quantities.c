/* MC-QUANTITIES.C - Module for quantities relating to Monte Carlo simulation.*/

/* Copyright (c) 1995, 1996 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 *
 * The 'k' quantity was implemented by Carl Edward Rasmussen, 1995.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "misc.h"
#include "log.h"
#include "quantities.h"
#include "mc.h"


/* STATE STRUCTURE. */

static mc_dynamic_state ds;


/* INITIALIZE AFTER FIRST RECORDS READ. */

void mc_initialize
( log_gobbled *logg
)
{ 
  mc_app_initialize(logg,&ds);

  if (logg->data['p'] && logg->actual_size['p'] != ds.dim * sizeof (mc_value))
  { fprintf(stderr,"Momentum record has wrong size (in mc_initialize.c)\n");
    exit(1);
  }

  logg->req_size['p'] = ds.dim * sizeof (mc_value);
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE. */

void mc_available 
( quantities_described qd,
  log_gobbled *logg
)
{ 
  int mod, low;
  int i;

  for (i = 0; i<Max_quantities; i++)
  {
    if (qd[i].letter && qd[i].available==0)
    {
      low = qd[i].low;
      mod = qd[i].modifier;

      if (strchr("TEKHDdfFmrekiIjJ",qd[i].letter)!=0)
      { qd[i].available = low==-1 && (mod==-1 || qd[i].letter=='D') ? 1 : -1;
      }

      if (qd[i].letter=='q')
      { qd[i].available = low==-1 && mod>=0 && mod<ds.dim ? 1 : -1;
      }

      if (qd[i].letter=='s')
      { qd[i].available = low==-1 && mod>=0 && mod<ds.dim ? 1 : -1;
      }
  
      if (qd[i].letter=='p')
      { qd[i].available = low==-1 && mod>=0 && mod<ds.dim ? 1 : -1;
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE. */

void mc_evaluate 
( quantities_described qd, 
  quantities_held *qh,
  log_gobbled *logg
)
{ 
  mc_iter *it;
  int mod;
  int i;
  
  it = logg->data['i'];
  ds.p = logg->data['p'];

  ds.know_grad    = 0;
  ds.know_pot     = 0;
  ds.know_kinetic = 0;
  
  for (i = 0; i<Max_quantities; i++)
  {
    if (qd[i].letter && !qh->updated[i])
    {
      mod = qd[i].modifier;

      switch (qd[i].letter)
      {
        case 'T':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = it->temperature;
          qh->updated[i] = 1;
          break;
        }

        case 'E':
        { if (!ds.know_pot)
          { mc_app_energy (&ds, 1, 1, &ds.pot_energy, (mc_value *) 0);
            ds.know_pot = 1;
          }
          *qh->value[i] = ds.pot_energy;
          qh->updated[i] = 1;
          break;
        }
  
        case 'K':
        { if (!ds.know_kinetic)
          { ds.kinetic_energy = mc_kinetic_energy(&ds);
            ds.know_kinetic = 1;
          }
          *qh->value[i] = ds.kinetic_energy;
          qh->updated[i] = 1;
          break;
        }
  
        case 'H':
        { if (!ds.know_pot)
          { mc_app_energy (&ds, 1, 1, &ds.pot_energy, (mc_value *) 0);
            ds.know_pot = 1;
          }
          if (!ds.know_kinetic)
          { ds.kinetic_energy = mc_kinetic_energy(&ds);
            ds.know_kinetic = 1;
          }
          *qh->value[i] = ds.pot_energy + ds.kinetic_energy;
          qh->updated[i] = 1;
          break;
        }
  
        case 'D':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = 
               it->delta > (mod<0?1000:mod) ? (mod<0?1000:mod) : it->delta;
          qh->updated[i] = 1;
          break;
        }

        case 'i': 
        { mc_temp_state *ts;
          ts = logg->data['b'];
          *qh->value[i] = ts==0 || logg->index['b']!=logg->last_index ? 1.0
                          : ts->inv_temp;
          qh->updated[i] = 1;
          break;
        }

        case 'I': 
        { 
          float inv_temp;
          mc_temp_sched *sch;
          mc_temp_state *ts;
          int j;

          if (logg->data['m']==0) break;

          ts = logg->data['b'];
          inv_temp = ts==0 || logg->index['b']!=logg->last_index ? 1.0
                     : ts->inv_temp;

          sch = logg->data['m'];
          for (j = 0; j<Max_temps && sch->sched[j].inv_temp!=1 
                       && sch->sched[j].inv_temp!=inv_temp; j++) ;
          if (sch->sched[j].inv_temp!=inv_temp) break; 

          *qh->value[i] = j;
          qh->updated[i] = 1;
          break;
        }

        case 'j': 
        { mc_temp_state *ts;
          ts = logg->data['b'];
          if (ts==0 || logg->index['b']!=logg->last_index) break; 
          *qh->value[i] = ts->temp_dir;
          qh->updated[i] = 1;
          break;
        }

        case 'J': 
        { float inv_temp;
          mc_temp_sched *sch;
          mc_temp_state *ts;
          int j, k;

          ts = logg->data['b'];
          if (ts==0 || logg->index['b']!=logg->last_index) break; 

          inv_temp = ts==0 || logg->index['b']!=logg->last_index ? 1.0
                     : ts->inv_temp;

          sch = logg->data['m'];

          for (j = 0; j<Max_temps && sch->sched[j].inv_temp!=1 
                       && sch->sched[j].inv_temp!=inv_temp; j++) ;
          if (sch->sched[j].inv_temp!=inv_temp) break; 

          k = j+ts->temp_dir;

          *qh->value[i] = j>k ? j : k;
          qh->updated[i] = 1;

          break;
        }
  
        case 'd':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = it->decay;
          qh->updated[i] = 1;
          break;
        }
  
        case 'f':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = it->stepsize_factor;
          qh->updated[i] = 1;
          break;
        }
  
        case 'F':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = it->stepsize_factor/2;
          qh->updated[i] = 1;
          break;
        }
  
        case 'm':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = it->move_point;
          qh->updated[i] = 1;
          break;
        }
  
        case 'r':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = (double)it->rejects / it->proposals;
          qh->updated[i] = 1;
          break;
        }

        case 'e':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = (double)it->slice_evals / it->slice_calls;
          qh->updated[i] = 1;
          break;
        }
  
        case 'q':
        { *qh->value[i] = ds.q[mod];
          qh->updated[i] = 1;
          break;
        }
    
        case 'p':
        { *qh->value[i] = ds.p==0 ? 0 : ds.p[mod];
          qh->updated[i] = 1;
          break;
        }
  
        case 's':
        { mc_app_stepsizes(&ds);
          *qh->value[i] = ds.stepsize[mod];
          qh->updated[i] = 1;
          break;
        }

        case 'k':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = (double)it->time/60000;
          qh->updated[i] = 1;
          break;
	}
      }
    }
  }
}





