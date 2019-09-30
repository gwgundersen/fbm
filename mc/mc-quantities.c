/* MC-QUANTITIES.C - Module for quantities relating to Monte Carlo simulation.*/

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
  char letter;
  int i;

  for (i = 0; i<Max_quantities; i++)
  {
    letter = qd[i].letter;
    low = qd[i].low;
    mod = qd[i].modifier;

    if (letter && qd[i].available==0 && (strchr("iIT",letter)==0 || mod==-1))
    {
      if (strchr("QEKHDdfmrekjJ",letter)!=0)
      { qd[i].available = 
          low==-1 && (mod==-1 || letter=='E' && (mod==0 || mod==1)
                              || letter=='K' && mod==0
                              || letter=='m' && mod<=5
#if 0
                              || letter=='H' && mod==0
#endif
                              || letter=='D' || letter=='Q')
          ? 1 : -1;
      }

      if (strchr("iIT",letter)!=0 && mod==-1)
      { qd[i].available = low==-1 ? 1 : -1;
      }

      if (letter=='F')
      { qd[i].available = low==-1 && (mod==1 || mod==2) ? 1 : -1;
      }

      if (letter=='q')
      { qd[i].available = low==-1 && mod>=0 && mod<ds.dim ? 1 : -1;
      }

      if (letter=='s')
      { qd[i].available = low==-1 && mod>=0 && mod<ds.dim ? 1 : -1;
      }
  
      if (letter=='p')
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
  char letter;
  int mod;
  int i;
  
  it = logg->data['i'];
  ds.p = logg->data['p'];

  ds.know_grad    = 0;
  ds.know_pot     = 0;
  ds.know_kinetic = 0;
  
  for (i = 0; i<Max_quantities; i++)
  {
    letter = qd[i].letter;
    mod = qd[i].modifier;

    if (letter && !qh->updated[i] && (strchr("iIT",letter)==0 || mod==-1))
    {
      switch (letter)
      {
        case 'T':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          *qh->value[i] = it->temperature;
          qh->updated[i] = 1;
          break;
        }

        case 'E':
        { if (!ds.know_pot && mod!=0)
          { mc_app_energy (&ds, 1, 1, &ds.pot_energy, (mc_value *) 0);
            ds.know_pot = 1;
          }
          if (mod==0 || mod==1)
          { mc_temp_state ts0, *ts1;
            ts0.inv_temp = 0;
            ts1 = ds.temp_state;
            ds.temp_state = &ts0;
            mc_app_energy (&ds, 1, 1, qh->value[i], 0);
            ds.temp_state = ts1;
          }
          if (mod==1)
          { *qh->value[i] = ds.pot_energy - *qh->value[i];
          }
          if (mod==-1)
          { *qh->value[i] = ds.pot_energy;
          }
          qh->updated[i] = 1;
          break;
        }
  
        case 'K':
        { if (mod==0)
          { *qh->value[i] = ds.dim/2.0;
            qh->updated[i] = 1;
          }
          else
          { if (!ds.know_kinetic)
            { ds.kinetic_energy = mc_kinetic_energy(&ds);
              ds.know_kinetic = 1;
            }
            *qh->value[i] = ds.kinetic_energy;
            qh->updated[i] = 1;
          }
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
#if 0
          if (mod==0)
          { mc_therm_state *h;
            h = logg->data['h'];
            if (h!=0 && logg->index['h']==logg->last_index)
            { *qh->value[i] += 0.5 * h->tp * h->tp;
              *qh->value[i] += log (h->tq>0 ? h->tq : -h->tq);
            }
          }
#endif 
          qh->updated[i] = 1;
          break;
        }

        case 'Q':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          if (mod<0)
          { *qh->value[i] = exp(it->log_weight);
          }
          else
          { *qh->value[i] = 
              mod==0 || it->log_weight>=-mod ? it->log_weight : -mod;
          }
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
          *qh->value[i] = mod==1 ? exp(it->log_tt_weight)
                                 : exp(it->log_tt_weight2);
          qh->updated[i] = 1;
          break;
        }
  
        case 'm':
        { if (it==0 || logg->index['i']!=logg->last_index) break;
          switch (mod)
          { case -1:
            { *qh->value[i] = it->move_point;
              break;
            }
            case 0:
            { *qh->value[i] = it->move_point + it->spiral_offset;
              break;
            }
            case 1:
            { *qh->value[i] = it->spiral_offset;
              break;
            }
            case 2:
            { *qh->value[i] = it->spiral_switch;
              break;
            }
            case 3:
            { *qh->value[i] = it->move_point + it->spiral_offset 
                                      - it->spiral_switch;
              break;
            }
            case 4:
            { *qh->value[i] = it->spiral_offset - it->spiral_switch;
              break;
            }
            case 5:
            { *qh->value[i] = 
                 abs(it->move_point + it->spiral_offset - it->spiral_switch)
               - abs(it->spiral_offset - it->spiral_switch);
              break;
            }
            default: abort();
          }
          qh->updated[i] = 1;
          break;
#if 0
          { mc_therm_state *h;
            h = logg->data['h'];
            if (h==0 || logg->index['h']!=logg->last_index) break;
            *qh->value[i] = mod==0 ? h->tq : h->tp;
            qh->updated[i] = 1;
          }
#endif
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





