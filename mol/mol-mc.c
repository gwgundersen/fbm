/* MOL-MC.C - Molecular dynamics module for Markov chain Monte Carlo. */

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


/* LOCAL VARIABLES. */

static mol_spec *ms;	/* Molecular dynamics specification */


/* COMPUTE ENERGY AND ITS GRADIENT FOR ONE PAIR OF PARTICLES.  The contribution
   from this pair is added to the variables pointed to by 'energy' and 'grad',
   if they are not null. */

static void ener_grad
( double *q,		/* Coordinates of molecules */
  double inv_temp,	/* Inverse temperature */
  mol_spec *ms,		/* Specifications of molecular system */
  int i,		/* Indexes of particles times */
  int j,
  double *energy,	/* Place to increment energy, or null */
  double *grad		/* Place to increment gradient, or null */
)
{
  double E, gri[3], grj[3], grl;
  int k, iD, jD;

  if (inv_temp==0) return;

  if (grad)
  { E = energy_and_gradient (q, inv_temp, ms, i, j, gri, grj, &grl);
  }
  else
  { E = energy_and_gradient (q, inv_temp, ms, i, j, 0, 0, 0);
  }

  if (energy) *energy += E;

  if (grad)
  { iD = i*ms->D;
    jD = j*ms->D;
    for (k = 0; k<ms->D; k++)
    { grad[iD+k] += gri[k];
      grad[jD+k] += grj[k];
    }
    if (ms->len_pres<0)
    { grad[ms->N*ms->D] += grl;
    }
  }
}


/* SET UP REQUIRED RECORD SIZES PRIOR TO GOBBLING RECORDS. */

void mc_app_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{ logg->req_size['M'] = sizeof *ms;
}


/* INITIALIZE AND SET UP DYNAMIC STATE STRUCTURE. */

void mc_app_initialize
( log_gobbled *logg,	/* Records gobbled up from head and tail of log file */
  mc_dynamic_state *ds	/* Structure holding pointers to dynamical state */
)
{ 
  double len;
  int i;

  if ((ms = logg->data['M'])==0)
  { fprintf(stderr,"No molecular dynamics specification in log file\n");
    exit(1);
  }

  ds->dim = ms->D*ms->N + (ms->len_pres<0);

  ds->temp_state = 0;

  logg->req_size['C'] = ds->dim * sizeof(mc_value);

  if ((ds->q = logg->data['C'])==0)
  { 
    ds->q = chk_alloc (ds->dim, sizeof(mc_value));
    if (ms->len_pres<0)
    { ds->q[ms->N*ms->D] = log(ms->width) + log(ms->N/0.5)/ms->D;
    }
    for (i = 0; i<ms->N*ms->D; i++) 
    { ds->q[i] = rand_uniform();
    }
  }
  else
  {
    if (logg->index['C'] != logg->last_index)
    { fprintf(stderr,"Missing current point in iteration stored in log file\n");
      exit(1);
    }
  }

  ds->stepsize = chk_alloc (ds->dim, sizeof(mc_value));
  mc_app_stepsizes(ds);
}


/* SAVE POSITION AND AUXILIARY PART OF STATE. */

void mc_app_save
( mc_dynamic_state *ds,	/* Current dyanamical state */
  log_file *logf,	/* Log file state structure */
  int index		/* Index of iteration being saved */
)
{
  logf->header.type = 'C';
  logf->header.index = index;
  logf->header.size = ds->dim * sizeof(mc_value);
  log_file_append(logf,ds->q);
}


/* APPLICATION-SPECIFIC SAMPLING PROCEDURES. */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double a,
  double a2,
  mc_iter *it,
  mc_temp_sched *sch
)
{
  if (strcmp(op,"wrap")==0)
  { int i;
    for (i = 0; i<ms->N*ms->D; i++)
    { ds->q[i] = wrap(ds->q[i]);
    }
    return 1;
  }

  if (strcmp(op,"met-mol")==0)
  { 
    double save[3];
    double Eold, Enew;
    double inv_temp;
    double U;
    int i, j, k;

    inv_temp = !ds->temp_state ? 1 : ds->temp_state->inv_temp;
    if (inv_temp<0) 
    { inv_temp = -inv_temp;
    }

    mc_app_stepsizes(ds);
    it->stepsize_factor = a;

    for (i = 0; i<ms->N; i++)
    { 
      Eold = 0;
      if (ms->scale!=0)
      { for (j = 0; j<ms->N; j++)
        { if (j!=i) 
          { ener_grad (ds->q,inv_temp,ms,i,j,&Eold,0);
          }
        }
      }

      for (k = 0; k<ms->D; k++) 
      { save[k] = ds->q[i*ms->D+k];
        ds->q[i*ms->D+k] += rand_gaussian() * ds->stepsize[0] * a;
      }

      Enew = 0;
      if (ms->scale!=0)
      { for (j = 0; j<ms->N; j++)
        { if (j!=i) 
          { ener_grad (ds->q,inv_temp,ms,i,j,&Enew,0);
          }
        }
      }

      it->proposals += 1;
      it->delta = Enew - Eold;

      U = rand_uniform();
      if (U<exp(-it->delta/it->temperature))
      { it->move_point = 1;
      }
      else
      { it->rejects += 1;
        it->move_point = 0;
        for (k = 0; k<ms->D; k++) 
        { ds->q[i*ms->D+k] = save[k];
        }
      }
    }

    ds->know_grad = ds->know_pot = 0;

    return 1;
  }

  if (strcmp(op,"met-len")==0)
  { 
    double Eold, Enew;
    double inv_temp;
    double s, U;
    int i, j, k;

    if (ms->len_pres>0)
    { fprintf(stderr,"The met-vol operation applies only to NPT ensembles\n");
      exit(1);
    }

    inv_temp = !ds->temp_state ? 1 : ds->temp_state->inv_temp;
    if (inv_temp<0) 
    { inv_temp = -inv_temp;
    }

    mc_app_stepsizes(ds);
    it->stepsize_factor = a;

    s = ds->q[ms->N*ms->D];

    mc_app_energy (ds, 1, 1, &Eold, 0);

    ds->q[ms->N*ms->D] = s + rand_gaussian() * ds->stepsize[ms->N*ms->D] * a;

    mc_app_energy (ds, 1, 1, &Enew, 0);

    it->proposals += 1;
    it->delta = Enew - Eold;

    U = rand_uniform();
    if (U<exp(-it->delta/it->temperature))
    { it->move_point = 1;
    }
    else
    { it->rejects += 1;
      it->move_point = 0;
      ds->q[ms->N*ms->D] = s;
    }

    ds->know_grad = ds->know_pot = 0;

    return 1;
  }

  return 0;
}


/* EVALUATE POTENTIAL ENERGY AND ITS GRADIENT. */

void mc_app_energy
( mc_dynamic_state *ds,	/* Current dyanamical state */
  int N_approx,		/* Number of gradient approximations in use */
  int w_approx,		/* Which approximation to use this time */
  double *energy,	/* Place to store energy, null if not required */
  mc_value *grad	/* Place to store gradient, null if not required */
)
{
  double inv_temp;
  int i, j, k, ninv;
  double P, L, E;

  inv_temp = !ds->temp_state ? 1 : ds->temp_state->inv_temp;
  ninv = 0;
  if (inv_temp<0) 
  { inv_temp = -inv_temp;
    ninv = 1;
  }

  /* Initialize energy and its gradient to zero. */

  if (grad) 
  { for (k = 0; k<ds->dim; k++) grad[k] = 0;
  }
 
  if (energy)
  { *energy = 0;
  }

  /* Look at each pair of particles, adding to the energy and gradient. */

  if (ms->scale!=0 && inv_temp!=0)
  {
    for (i = 0; i<ms->N; i++)
    { for (j = i+1; j<ms->N; j++)
      { ener_grad (ds->q, inv_temp, ms, i, j, energy, grad);
      }
    }
  }

  /* Consider potential energy term relating to volume, for NPT ensemble. 
     If the log length of a dimension is more than 100 or less than -100, 
     consider the energy to be as if it was at this extreme, and let the
     gradient be zero. */

  if (ms->len_pres<0)
  {
    P = -ms->len_pres;
    L = ds->q[ms->N*ms->D];
    if (L>100)
    { L = 100;
    }
    if (L<-100)
    { L = -100;
    }

    if (ninv)
    { inv_temp = 1-ms->inv_temp;
    }
    else if (inv_temp==0)
    { inv_temp = ms->inv_temp;
    }

    if (energy)
    { *energy += inv_temp * P*exp(L*ms->D);
      if (!ninv)
      { *energy -= (ms->N+1)*L*ms->D;
      }
    }
    if (grad && ds->q[ms->N*ms->D]<100 && ds->q[ms->N*ms->D]>-100)
    { grad[ms->N*ms->D] += inv_temp * ms->D*P*exp(L*ms->D) - (ms->N+1)*ms->D;
      if (!ninv)
      { grad[ms->N*ms->D] -= (ms->N+1)*ms->D;
      }
    }
  }

  if (energy && isnan(*energy)) abort();
  if (grad)
  { for (k = 0; k<ds->dim; k++)
    { if (isnan(grad[k])) abort();
    }
  }
}


/* SAMPLE FROM DISTRIBUTION AT INVERSE TEMPERATURE OF ZERO.  Returns zero
   if this is not possible.  For NPT ensemble, the "infinite temperature" 
   distribution for volume isn't really infinite temperature, but has 
   the temperature specified in mol-spec. */

int mc_app_zero_gen
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
  double vol, len;
  int i;

  if (ms->len_pres<0)
  { vol = rand_gamma(ms->N+1.0) / (-ms->len_pres*ms->inv_temp);
    ds->q[ms->N*ms->D] = log(vol) / ms->D;
  }

  for (i = 0; i<ms->N*ms->D; i++) 
  { ds->q[i] = rand_uniform();
  }

  return 1;
}


/* SET STEPSIZES FOR EACH COORDINATE. */

void mc_app_stepsizes
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
  double f, s;
  int i;

  f = ds->temp_state==0 || ds->temp_state->inv_temp==0 ? 1 
        : 1 / sqrt(fabs(ds->temp_state->inv_temp));

  s = f * pow(ms->N,-1.0/ms->D);

  for (i = 0; i<ms->N*ms->D; i++) 
  { ds->stepsize[i] = s;
  }

  if (ms->len_pres<0)
  { ds->stepsize[ms->N*ms->D] = f / ms->N;
  }
}
