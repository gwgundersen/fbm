/* MC-ITER.C - Procedures for performing Markov chain Monte Carlo iterations. */

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
#include "quantities.h"


/* OPTIONS. */

#define Print_proposals 0	/* Print proposed changes by quadrant? 
				   (not normally done, for test only) */

static int echo = 0;		/* Debug option - echos operations first time */

#define Max_temp_repeat	1000	/* Max repeat count in a tempered transition */


/* LOCAL VARIABLES. */

static mc_ops *ops;	/* Operations to perform each iteration */
static int endo;               /* Index of point past last operation */
static int endp[2*Max_mc_ops]; /* Other endpoint for group operations */

static mc_traj *tj;	/* Trajectory specification */
static mc_temp_sched *sch; /* Tempering schedule, zero if none */

static int need_p;	/* Do we need momentum? */
static int need_grad;	/* Do we need gradient? */

static int have_ss;	/* Do we have the stepsizes? */

static int need_save;	/* Do we need space to save q and (maybe) p? */
static mc_value *q_save;/* Place to save q */
static mc_value *p_save;/* Place to save p */

static int need_lowhigh;/* Do we need space for low and high bounds? */
static mc_value *lowb;	/* Place to store low bound (for "slice") */
static mc_value *highb;	/* Place to store high bound (for "slice") */

static int need_wsum;	/* Do we need space to save weighted sum of crumbs? */
static mc_value *wsum;	/* Place to store wsum (for "slice-gaussian") */

static int need_savet;	/* Do we need space to save for tempered transitions? */
static mc_value *q_savet;/* Place to save q for tempered transitions */
static mc_value *p_savet;/* Place to save p for tempered transitions */
static mc_value *aux_savet; /* Place to save auxiliary variables */

static int need_arsv;	/* Do we need space to save accept/reject states? */
static mc_value *q_asv;	/* Place to save q values for accept state */
static mc_value *p_asv;	/* Place to save p values for accept state */
static mc_value *q_rsv;	/* Place to save q values for reject state */
static mc_value *p_rsv;	/* Place to save p values for reject state */

static int print_index;	/* Index used to label printed quantities */


/* LOCAL PROCEDURES.  Tempering procedures are included in this module
   to simplify the interface. */

static void do_group (mc_dynamic_state *, mc_iter *, int, int, int,
                      log_gobbled *, quantities_described, int);

void mc_simulated_tempering (mc_dynamic_state *, mc_iter *);
void mc_tempered_transition (mc_dynamic_state *, mc_iter *, 
                             int, int, int, log_gobbled *, 
                             quantities_described, int);


/* STRUCTURE AND PROCEDURE USED FOR SORTING. */

typedef struct { int index; double value; } sort_table;

int compare (const void *a, const void *b) 
{ 
  return ((sort_table*)a)->value < ((sort_table*)b)->value ? -1
       : ((sort_table*)a)->value > ((sort_table*)b)->value ? +1
       : 0;
}


/* SET UP FOR PERFORMING ITERATIONS. */

void mc_iter_init 
( mc_dynamic_state *ds,	/* State to update */
  mc_ops *ops0,		/* Operations to perform each iteration */
  mc_traj *tj0,		/* Trajectory specification */
  mc_temp_sched *sch0	/* Tempering schedule, zero if none */
)
{  
  int depth, type, i;
  int stack[Max_mc_ops];
  int does_print;

  ops = ops0;
  tj = tj0;
  sch = sch0;

  need_p = need_grad = need_save = need_lowhigh = need_wsum =
           need_savet = need_arsv = 0;

  does_print = 0;

  depth = 0;

  for (i = 0; i<Max_mc_ops && ops->op[i].type!=0; i++)
  { 
    type = ops->op[i].type;

    endp[i] = i;

    if (strchr(Group_ops,type)!=0)
    { stack[depth] = i;
      depth += 1;
    }

    if (type=='E')
    { depth -= 1;
      if (depth<0) 
      { fprintf(stderr,"Too many 'end' ops in Monte Carlo operations list\n");
        exit(1);
      }
      endp[stack[depth]] = i;
      endp[i] = stack[depth];
    }

    if (strchr(MC_needs_p,type)!=0)
    { need_p = 1;
    }

    if (strchr(MC_needs_grad,type)!=0 && (type!='l' || ops->op[i].g_shrink!=0))
    { need_grad = 1;
    }

    if (type=='M' || type=='H' || type=='T' || type=='@' || type=='^'
     || type=='i' || type=='o' || type=='G' || type=='l' || type=='u')
    { need_save = 1;
    }

    if (type=='H' || type=='T' || type=='@' || type=='^')
    { need_arsv = 1;
    }

    if (type=='t')
    { need_savet = 1;
    }

    if (type=='l')
    { need_lowhigh = 1;
    }

    if (type=='u')
    { need_wsum = 1;
    }

    if (type=='p') 
    { does_print = 1;   
    }
  }

  while (depth>0)
  { depth -= 1;
    endp[stack[depth]] = i;
    endp[i] = stack[depth];
    i += 1;
  }

  endo = i-1;

  if (need_save)
  { q_save = chk_alloc (ds->dim, sizeof *q_save);
    p_save = chk_alloc (ds->dim, sizeof *p_save);
  }

  if (need_lowhigh)
  { lowb  = chk_alloc (ds->dim, sizeof *lowb);
    highb = chk_alloc (ds->dim, sizeof *highb);
  }

  if (need_wsum)
  { wsum = chk_alloc (ds->dim, sizeof *lowb);
  }

  if (need_savet)
  { q_savet = chk_alloc (ds->dim, sizeof *q_savet);
    p_savet = chk_alloc (ds->dim, sizeof *p_savet);
    if (ds->aux_dim<0)
    { fprintf(stderr,
        "Saving of auxiliary state isn't implemented (temp-trans won't work)\n");
      exit(1);
    }
    if (ds->aux_dim>0) aux_savet = chk_alloc (ds->aux_dim, sizeof *p_savet);
  }

  if (need_arsv)
  { q_asv = chk_alloc (ds->dim, sizeof *q_asv);
    p_asv = chk_alloc (ds->dim, sizeof *p_asv);
    q_rsv = chk_alloc (ds->dim, sizeof *q_rsv);
    p_rsv = chk_alloc (ds->dim, sizeof *p_rsv);
  }
}


/* PERFORM ONE ITERATION.  A pointer to the structure holding data on this
   iteration is passed.  The caller will have set the temperature and decay 
   fields as specified by the user.  The other fields are set by this 
   procedure, or added to, in the case of rejects and proposals. */

void mc_iteration
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  log_gobbled *logg,	/* Records gobbled */
  void *qd, 		/* Descriptions of quantities to plot */
  int N_quantities	/* Number of quantities to plot, -1 for tt plot */
)
{ 
  int j, na;

  /* Create momentum variables and gradient vector if needed. */

  if (need_p && ds->p==0)
  { 
    ds->p = chk_alloc (ds->dim, sizeof (mc_value));
    mc_heatbath (ds, it->temperature, 0.0);
    ds->know_kinetic = 0;
  }

  if (need_grad && ds->grad==0)
  {
    ds->grad = chk_alloc (ds->dim, sizeof (mc_value));
    ds->know_grad = 0;
  }

  /* Initialize fields describing iteration, except those that are additive. */

  it->stepsize_factor = 1.0;
  it->move_point = 0;
  it->delta = 0.0;
  it->log_tt_weight = 1e30;
  it->log_tt_weight2 = 1e30;

  /* Set approximation order, usually superseded by a permute call. */

  na = tj->N_approx>0 ? tj->N_approx : -tj->N_approx;
  if (na>Max_approx) na = Max_approx;

  for (j = 0; j<na; j++) it->approx_order[j] = j+1;

  /* Perform the operations. */

  have_ss = 0;

  if (echo)
  { printf("Ops performed:");
  }
 
  print_index = 0;

  do_group(ds,it,0,endo,0,logg,qd,N_quantities);

  if (echo)
  { printf("\n");
    echo = 0;  
  }
}


/* PERFORM A GROUP OF OPERATIONS. */

static void do_group
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int start,		/* Index of first operation in group */
  int end,		/* Index of terminator for group of operations */
  int reverse,		/* Do operations in reverse order? */
  log_gobbled *logg,	/* Records gobbled */
  quantities_described qd, /* Descriptions of quantities to plot */
  int N_quantities	/* Number of quantities to plot, -1 for tt plot */
)
{
  double alpha, stepsize_adjust;
  int type;

  int i, k, c;

  if (end<start) return;

  if (reverse)
  { i = endp[end];
  }
  else
  { i = start;
  }

  for (;;)
  {
    if (i>=Max_mc_ops) abort();

    type = ops->op[i].type;

    if (type==0) abort();

    if (echo)
    { printf(" %d%c",i,type);
    }

    /* Figure out what stepsize factor to use this time, if necessary. */

    if (type=='M' || type=='m' || type=='S' || type=='O' 
     || type=='D' || type=='P' || type=='H' || type=='T'
     || type=='@' || type=='^' || type=='i' || type=='o' 
     || type=='h' || type=='G' || type=='g' || type=='l'
     || type=='u')
    { 
      stepsize_adjust = ops->op[i].stepsize_adjust;

      alpha = ops->op[i].stepsize_alpha;

      it->stepsize_factor = 
        stepsize_adjust>0 ? stepsize_adjust : -stepsize_adjust;

      if (alpha!=0) 
      { it->stepsize_factor *= 
          alpha>0 ? 1 / sqrt (rand_gamma(alpha/2) / (alpha/2))
                  : pow(10.0,-alpha*(rand_uniopen()-0.5));
      }

      if (stepsize_adjust<0)
      { for (k = 0; k<ds->dim; k++) ds->stepsize[k] = 1.0;
        have_ss = 0;
      }
      else if (!have_ss)
      { mc_app_stepsizes (ds);
        have_ss = 1;
      }
    }

    /* Do the next operation. */

    switch (type)
    { 
      case 'E': abort();

      case 'R': 
      { for (c = 0; c<ops->op[i].repeat_count; c++)
        { do_group(ds,it,i+1,endp[i]-1,reverse,logg,qd,N_quantities);
        }
        break;
      }

      case 'x':
      { if (!have_ss)
        { mc_app_stepsizes (ds);
          have_ss = 1;
        }
        for (k = (ops->op[i].firsti==-1 ? 0 : ops->op[i].firsti); 
             k <= (ops->op[i].firsti==-1 ? ds->dim : ops->op[i].lasti)
              && k<ds->dim;
             k++)
        { ds->stepsize[k] *= ops->op[i].stepsize_adjust;
        }
        break;
      }

      case 'B':
      { double d;
        d = it->decay>=0 ? it->decay : ops->op[i].heatbath_decay;
        mc_heatbath (ds, it->temperature, d);
        break;
      }

      case 'r':
      { mc_radial_heatbath (ds, it->temperature);
        break;
      }

      case 'X':
      { mc_mix_momentum (ds, ops->op[i].heatbath_decay);
        break;
      }

      case 'N':
      { for (k = 0; k<ds->dim; k++) 
        { ds->p[k] = -ds->p[k];
        }
        break;
      }

      case 'M': 
      { mc_metropolis(ds,it,q_save,ops->op[i].b_accept);
        break;
      }

      case 'G': 
      { mc_rgrid_met(ds,it,q_save,ops->op[i].b_accept);
        break;
      }

      case 'm':
      { mc_met_1(ds,it,ops->op[i].firsti,ops->op[i].lasti,
                 ops->op[i].b_accept, ops->op[i].r_update);
        break;
      }

      case 'g':
      { mc_rgrid_met_1(ds,it,ops->op[i].firsti,ops->op[i].lasti,
                       ops->op[i].b_accept, ops->op[i].r_update);
        break;
      }

      case 'U':
      { mc_gaussian_gibbs(ds,it,ops->op[i].firsti,ops->op[i].lasti,
                          ops->op[i].r_update);
        break;
      }

      case 'D':
      { mc_traj_init(tj,it);
        mc_trajectory(ds,ops->op[i].steps,0);
        break;
      }
#if 0
      case 'h':
      { mc_therm_present(ds);
        mc_traj_init(tj,it);
        mc_therm_trajectory(ds,ops->op[i].steps,0);
        break;
      }
#endif
      case 'P':
      { mc_traj_init(tj,it);
        mc_traj_permute();
        mc_trajectory(ds,ops->op[i].steps,0);
        break;
      }

      case 'H': case 'T':
      { if (ops->op[i].in_steps==0)
        { mc_hybrid (ds, it, tj, ops->op[i].steps, ops->op[i].window, 
                  ops->op[i].jump, type=='H' ? 0.0 : ops->op[i].temper_factor, 
                  N_quantities, q_save, p_save, q_asv, p_asv, q_rsv, p_rsv);
        }
        else
        { mc_hybrid2 (ds, it, tj, ops->op[i].steps, ops->op[i].in_steps,
                ops->op[i].jump, q_save, p_save);
        }
        break;
      }

      case '@': case '^':
      { mc_spiral (ds, it, tj, ops->op[i].steps, ops->op[i].temper_factor, 
                   type=='^', q_save, p_save, q_asv, p_asv);
        break;
      }

      case 'S':
      { mc_slice_1 (ds, it, ops->op[i].firsti, ops->op[i].lasti, 
                    ops->op[i].steps, ops->op[i].r_update, ops->op[i].s_factor,
                    ops->op[i].s_threshold);
        break;
      }

      case 'l':
      { mc_slice (ds, it, q_save, lowb, highb, ops->op[i].g_shrink);
        break;
      }

      case 'u':
      { mc_slice_gaussian (ds, it, q_save, wsum, ops->op[i].e_shrink);
        break;
      }

      case 'O':
      { mc_slice_over (ds, it, ops->op[i].refinements, ops->op[i].refresh_prob,
                       ops->op[i].firsti, ops->op[i].lasti, ops->op[i].steps,
                       ops->op[i].r_update);
        break;
      }

      case 'i':
      { mc_slice_inside (ds, it, ops->op[i].steps, q_save, p_save);
        break;
      }
 
      case 'o':
      { mc_slice_outside (ds, it, ops->op[i].steps, ops->op[i].in_steps, 
                          q_save, p_save);
        break;
      }

      case 's':
      { mc_simulated_tempering(ds,it);
        have_ss = 0;
        break;
      }

      case 't':
      { mc_tempered_transition (ds, it, i+1, endp[i]-1, reverse, 
          logg, qd, N_quantities);
        break;
      }

      case 'n':
      { mc_temp_present(ds,sch);
        ds->temp_state->temp_dir = -ds->temp_state->temp_dir;
        break;
      }

      case 'b':
      { mc_temp_present(ds,sch);
        ds->temp_state->temp_dir = rand_int(2) ? -1 : +1;
        break;
      }

      case '*':
      { int j;
        double f;
        f = ops->op[i].heatbath_decay+1;
        for (j = 0; j<ds->dim; j++)
        { ds->p[j] *= f;
        }
        ds->know_kinetic = 0;
        break;
      }

      case '=':
      { int j;
        double v;
        v = ops->op[i].heatbath_decay;
        for (j = 0; j<ds->dim; j++)
        { ds->p[j] = v;
        }
        ds->know_kinetic = 0;
        break;
      }

      case 'a':
      {
        double old_energy;

        mc_temp_present(ds,sch);

        if (ds->temp_index<0 || ds->temp_index>Max_temps) abort();

        if (ds->temp_state->inv_temp==1)
        { if (!mc_app_zero_gen(ds))
          { fprintf(stderr,"Application doesn't support use of AIS\n");
            exit(1);
          }
          ds->temp_state->inv_temp = 0;
          mc_app_energy (ds, 1, 1, &old_energy, 0);
          ds->temp_index = 0;
          it->log_weight = 0;
        }
        else
        { if (!ds->know_pot) mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
          old_energy = ds->pot_energy;
          ds->temp_index += 1;
        }

        ds->temp_state->inv_temp = sch->sched[ds->temp_index].inv_temp;

        mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
        ds->know_pot = 1;
        ds->know_grad = 0;

        it->log_weight += old_energy - ds->pot_energy;

        have_ss = 0;

        break;
      }

      case 'p':
      { quantities_held *qh;
        int v;

        if (N_quantities<=0) break;

        if (N_quantities==1) 
        { if (print_index==0) printf("\n");
          printf("%6d",print_index);
        }

        qh = quantities_storage(qd);
        quantities_evaluate(qd,qh,logg);
        for (v = 0; v<N_quantities; v++) 
        { printf(" %20.8e",*qh->value[v]);
        }
        printf("\n");
        print_index += 1;
        free(qh);

        break;
      }

      case 'A':
      { if (!mc_app_sample (ds, ops->op[i].appl, ops->op[i].app_param, 
                            ops->op[i].app_param2, it, sch))
        { fprintf(stderr,"Unknown application-specific operation: %s\n",
            ops->op[i].appl);
          exit(1);
        }
        have_ss = 0;
        break;
      }

      default:
      { fprintf(stderr,"Unknown operation type encountered: %c\n", type);
        exit(1);
      }
    }

    if (reverse && i==start || !reverse && endp[i]==end) break;

    if (reverse)
    { i = endp[i-1];
    }
    else
    { i = endp[i]+1;
    }
      
  }
}


/* PERFORM SIMULATED TEMPERING UPDATE.  Does a Metropolis update for a
   proposal to increase or decrease the inverse temperature. */

void mc_simulated_tempering
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it		/* Description of this iteration */
)
{
  double old_energy, olde, newe, U;

  mc_temp_present(ds,sch);

  if (ds->temp_index<0 || ds->temp_index>Max_temps
   || ds->temp_state->temp_dir!=-1 && ds->temp_state->temp_dir!=+1)
  { abort();
  }

  if (!ds->know_pot)
  { mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);
    ds->know_pot = 1;
  }

  old_energy = ds->pot_energy;
  olde = old_energy + sch->sched[ds->temp_index].bias;

  if (ds->temp_index==0 && ds->temp_state->temp_dir==-1 
   || ds->temp_state->inv_temp==1 && ds->temp_state->temp_dir==+1)
  { it->proposals += 1;
    it->delta = 1e30;
    it->move_point = 0;
    it->rejects += 1;
    return;
  }

  ds->temp_index += ds->temp_state->temp_dir;
  ds->temp_state->inv_temp = sch->sched[ds->temp_index].inv_temp;
  
  mc_app_energy (ds, 1, 1, &ds->pot_energy, 0);

  newe = ds->pot_energy + sch->sched[ds->temp_index].bias;

  it->proposals += 1;
  it->delta = newe - olde;

  U = rand_uniform(); /* Do every time to keep in sync for coupling purposes */

  if (it->delta<=0 || U<exp(-it->delta/it->temperature))
  { 
    it->move_point = 1;
    ds->temp_state->temp_dir = -ds->temp_state->temp_dir;
    ds->know_grad = 0;
  }
  else
  { 
    it->move_point = 0;
    it->rejects += 1;
    ds->pot_energy = old_energy;
    ds->temp_index -= ds->temp_state->temp_dir;
    ds->temp_state->inv_temp = sch->sched[ds->temp_index].inv_temp;
  }
  
}


/* DO A TEMPERED TRANSITION. */

void mc_tempered_transition
( mc_dynamic_state *ds,	/* State to update */
  mc_iter *it,		/* Description of this iteration */
  int start,		/* Index of first operation in group */
  int end,		/* Index of terminator for group of operations */
  int reverse,		/* Do operations in reverse order? */
  log_gobbled *logg,	/* Records gobbled */
  quantities_described qd, /* Descriptions of quantities to plot */
  int N_quantities	/* Number of quantities to plot, -1 for tt plot */
)
{
  double down[Max_temps];
  mc_temp_state ts;
  double delta, ed, b1, b2, U;
  int quad0, quad1;
  int i1, i2;

  if (ds->temp_state)
  { fprintf(stderr, "Tempering operations can't be nested\n");
    exit(1);
  }
  ds->temp_state = &ts;

  if (sch==0)
  { fprintf(stderr,"No tempering schedule has been specified\n");
    exit(1);
  }

  if (Print_proposals)
  { quad0 = ds->q[0]>0 && ds->q[1]>0 ? 1
          : ds->q[0]<=0 && ds->q[1]>0 ? 2
          : ds->q[0]<=0 && ds->q[1]<=0 ? 3
          : 4;
  }

  if (logg->data['b']!=0) abort();
  logg->data['b'] = ds->temp_state;
  logg->index['b'] = logg->last_index;

  ds->temp_state->inv_temp = 1.0;
  ds->temp_index = mc_temp_index(sch,1.0);
  ds->temp_state->temp_dir = -1;

  if (N_quantities<0)
  { printf("\n");
  }

  mc_value_copy (q_savet, ds->q, ds->dim);
  if (ds->p) mc_value_copy(p_savet, ds->p, ds->dim);
  if (ds->aux) mc_value_copy(aux_savet, ds->aux, ds->aux_dim);

  delta = 0;

  for (;;)
  { 
    b1 = ds->temp_state->inv_temp;
    i1 = ds->temp_index;

    if (ds->temp_state->temp_dir==-1)
    {
      if (ds->temp_index==0)
      { ds->temp_state->temp_dir = +1;
        it->log_tt_weight = -delta;
      }
      else
      { ed = mc_energy_diff(ds,sch,ds->temp_state->temp_dir);
        delta += ed;
        ds->temp_index -= 1;
        ds->know_pot = 0;
        ds->know_grad = 0;
        have_ss = 0;
      }
    }
    else
    { 
      ed = mc_energy_diff(ds,sch,ds->temp_state->temp_dir);
      delta += ed;
      ds->temp_index += 1;
      ds->know_pot = 0;
      ds->know_grad = 0;
      have_ss = 0;
    }

    b2 = ds->temp_state->inv_temp = sch->sched[ds->temp_index].inv_temp;
    i2 = ds->temp_index;

    if (N_quantities==-1) 
    { if (b1==b2)
      { printf("\n");
      }
      else 
      { printf(" %12.10f %20.10e\n", b1, ed/(b2-b1));
        printf(" %12.10f %20.10e\n", b2, ed/(b2-b1));
      }
    }

    if (N_quantities==-2)
    { if (ds->temp_state->temp_dir==-1)
      { down[i1] = ed/(b2-b1);
      }
      else
      { printf(" %12.10f %20.10e\n", b1, ed/(b2-b1) - down[i2]);
        printf(" %12.10f %20.10e\n", b2, ed/(b2-b1) - down[i2]);
      }
    }
    
    if (ds->temp_state->inv_temp==1.0) 
    { break;  
    }

    do_group (ds, it, start, end, 
              (ds->temp_state->temp_dir==-1 ? reverse : !reverse), 
              logg, qd, (N_quantities<0 ? 0 : N_quantities));
  }

  if (Print_proposals)
  { quad1 = ds->q[0]>0 && ds->q[1]>0 ? 1
          : ds->q[0]<=0 && ds->q[1]>0 ? 2
          : ds->q[0]<=0 && ds->q[1]<=0 ? 3
          : 4;
  }

  it->proposals += 1;
  it->delta = delta;

  U = rand_uniform(); /* Do every time to keep in sync for coupling purposes */

  if (delta<=0 || U<exp(-delta/it->temperature))
  { 
    it->move_point = 1;

    it->log_tt_weight2 = 
      addlogs (it->log_tt_weight, delta+it->log_tt_weight) - log(2.0);
  }
  else
  { 
    it->rejects += 1;
    it->move_point = 0;

    mc_value_copy (ds->q, q_savet, ds->dim);
    if (ds->p) mc_value_copy(ds->p, p_savet, ds->dim);
    if (ds->aux) mc_value_copy(ds->aux, aux_savet, ds->aux_dim);

    it->log_tt_weight2 = it->log_tt_weight;
  }

  if (Print_proposals)
  { printf("%d %d %d %.1f\n",quad0,quad1,it->move_point,it->delta);
    fflush(stdout);
  }

  logg->data['b'] = 0;
  logg->index['b'] = 0;

  ds->temp_state = 0;
  ds->know_pot = 0;
  ds->know_grad = 0;
  have_ss = 0;
}
