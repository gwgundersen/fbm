/* DFT-MC.C - Interface between diffusion tree and Markov chain modules. */

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

#ifdef FLT_EXCEPTIONS
#include <fenv.h>
#endif

#include "misc.h"
#include "rand.h"
#include "log.h"
#include "mc.h"
#include "data.h"
#include "prior.h"
#include "model.h"
#include "dft.h"
#include "dft-data.h"


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* DIFFUSION TREE MODEL VARIABLES. */

static int initialize_done = 0;	/* Has this all been set up? */

static dft_spec *dft;		/* Diffusion tree model specification */
static model_specification *model; /* Data model */

static int N_targets;		/* Number of target values in dataset */

static dft_hypers *hyp;		/* Hyperparameters for diffusion tree model */

static int have_trees;		/* Are the trees present? */
static int have_latent;		/* Are latent vectors present? */
static int have_locations;	/* Are the node locations present? */

static int *parents;		/* Parents of nodes */
static double *divt;		/* Divergence times */
static double *latent;		/* Latent vectors */
static double *locations;	/* Node locations */

static dft_tree_node *nodes;	/* Internal tree representation */

static double *terminals[Max_trees]; /* Terminal node locations (temporary) */

static dft_state st;		/* Pointers to components of state */

static double *saved_location;	/* Place to temporarily save location of node */

static int *flags;	        /* Flags indexed by numbers from 0 to N_train */
static int flags_n;		/* No. of training cases flags allocated for */

static double inv_temp;		/* Inverse temperature */


/* PROCEDURES. */

static void create_trees (void), create_latent (void), create_locations (void), 
            gibbs_noise (void), gibbs_d_SD (void), 
            slice_div (mc_iter *, double, int),
            mh_latent (mc_iter *, double, int),
            sample_latent (void), sample_locations (void),
            slice_positions (mc_iter *, int, int), 
            met_terminals (mc_iter *, int, int),
            met_nonterminals (mc_iter *, int);


/* SET UP REQUIRED RECORD SIZES PRIOR TO GOBBLING RECORDS. */

void mc_app_record_sizes
( log_gobbled *logg	/* Structure to hold gobbled data */
)
{ 
  dft_record_sizes(logg);
}


/* INITIALIZE AND SET UP DYNAMIC STATE STRUCTURE.  Skips some stuff
   if it's already been done, as indicated by the initialize_done
   variable.  The real dynamic state structure is currently empty. */

void mc_app_initialize
( log_gobbled *logg,	/* Records gobbled up from head and tail of log file */
  mc_dynamic_state *ds	/* Structure holding pointers to dynamical state */
)
{ 
  int dt;

  if (!initialize_done)
  {
    /* Enable floating-point exceptions on invalid operations, if FLT_EXCEPTIONS
       is defined. */

#   ifdef FLT_EXCEPTIONS
      feenableexcept(FE_INVALID);
#   endif

    /* Check that required specification records are present. */
  
    dft   = logg->data['P'];
    model = logg->data['M'];

    dft_check_specs_present (dft, 0, model);

    N_targets = dft->N_targets;

    /* Check model type is OK. */

    if (model!=0 && model->type!='R' && model->type!='B')
    { fprintf(stderr,
      "The data model used is not implemented for Dirichlet diffusion trees\n");
      exit(1);
    }
  
    /* Locate existing state records. */
  
    hyp = logg->data['S'];

    have_trees     = logg->data['T']!=0;
    have_latent    = logg->data['L']!=0;
    have_locations = logg->data['N']!=0;

    if (have_trees && logg->data['R']==0)
    { fprintf(stderr,"Missing record of parents!\n");
      exit(1);
    }

    if (!have_trees && logg->data['R']!=0)
    { fprintf(stderr,"Missing record of divergence times!\n");
      exit(1);
    }
  
    if (hyp==0 && (have_trees || have_latent || have_locations))
    { fprintf(stderr,"Missing hyperparameter record!\n");
      exit(1);
    }

    if (!have_trees && (have_latent || have_locations))
    { fprintf(stderr,"Missing records describing tree structures!\n");
      exit(1);
    }

    if (hyp==0) 
    { hyp = chk_alloc (1, dft_hyper_size(dft,model));
      dft_hyper_init (dft, model, hyp);
    }

    /* Read training data, if any. */
  
    data_spec = logg->data['D'];

    if (data_spec && logg->actual_size['D'] !=
                       data_spec_size(data_spec->N_inputs,data_spec->N_targets))
    { fprintf(stderr,"Data specification record is the wrong size!\n");
      exit(1);
    }

    if (data_spec!=0) 
    { dft_data_read (1, 0, dft, model);
    }
    else
    { N_train = 0;
    }

    /* Check sizes of records holding state info. */

    dft_check_sizes(logg,dft,model,N_train);

    /* Set up pointers to existing state, or to newly allocated space. */

    parents   = have_trees    ? logg->data['R'] 
                              : chk_alloc(1,dft_parents_size(dft,N_train));
    divt      = have_trees    ? logg->data['T']
                              : chk_alloc(1,dft_divt_size(dft,N_train));
    latent    = have_latent   ? logg->data['L']
                              : chk_alloc(1,dft_latent_size(dft,N_train));
    locations = have_locations? logg->data['N']
                              : chk_alloc(1,dft_locations_size(dft,N_train));

    for (dt = 0; dt<dft->N_trees; dt++)
    { terminals[dt] = chk_alloc (N_train*N_targets, sizeof (double));
    }

    /* Make sure we don't do all this again. */

    initialize_done = 1;
  }

  /* Set up Monte Carlo state structure.  Everything is null. */

  ds->dim = 0;
  ds->q   = 0;

  ds->aux_dim = -1;  /* Indicates that tempered transitions work */
  ds->aux = 0;

  ds->temp_state = 0;
  
  ds->stepsize = 0;
}


/* RESET INITIALIZE_DONE IN PREPARATION FOR NEW LOG FILE. */

void dft_mc_cleanup(void)
{
  initialize_done = 0;
}


/* SAVE POSITION AND AUXILIARY PART OF STATE. */

void mc_app_save
( mc_dynamic_state *ds,	/* Current dyanamical state */
  log_file *logf,	/* Log file state structure */
  int index		/* Index of iteration being saved */
)
{ 
  logf->header.type = 'S';
  logf->header.index = index;
  logf->header.size = dft_hyper_size(dft,model);
  log_file_append (logf, hyp);

  if (N_train>0)
  {
    if (have_trees)
    { 
      logf->header.type = 'T';
      logf->header.index = index;
      logf->header.size = dft_divt_size(dft,N_train);
      log_file_append (logf, divt);

      logf->header.type = 'R';
      logf->header.index = index;
      logf->header.size = dft_parents_size(dft,N_train);
      log_file_append (logf, parents);
    }
  
    if (have_latent)
    { logf->header.type = 'L';
      logf->header.index = index;
      logf->header.size = dft_latent_size(dft,N_train);
      log_file_append (logf, latent);
    }
  
    if (have_locations)
    { logf->header.type = 'N';
      logf->header.index = index;
      logf->header.size = dft_locations_size(dft,N_train);
      log_file_append (logf, locations);
    }
  }
}


/* APPLICATION-SPECIFIC SAMPLING PROCEDURE. */

int mc_app_sample 
( mc_dynamic_state *ds,
  char *op,
  double pm,
  double pm2,
  mc_iter *it,
  mc_temp_sched *sch
)
{
  int k, K;

  inv_temp = !ds->temp_state ? 1 : ds->temp_state->inv_temp;

  /* Can't trust things when the inverse temperature is zero. */

  if (inv_temp==0) 
  { 
    fprintf(stderr,"Can't do operations at inverse temperature zero\n");
    exit(1);
  }

  ds->know_pot = ds->know_grad = 0;

  create_trees();

  dft_setup_state (dft, model, st, hyp, parents, divt, 
                   have_locations ? locations : 0, nodes, N_train);

  if (saved_location==0)
  { saved_location = chk_alloc (dft->N_targets, sizeof *saved_location);
  }

  if (strcmp(op,"gibbs-noise")==0 || strcmp(op,"gibbs-hypers")==0
   || strcmp(op,"gibbs-sigmas")==0)
  { 
    int had_locations, had_latent, do_noise, do_hypers;

    do_noise  = strcmp(op,"gibbs-hypers")!=0;
    do_hypers = strcmp(op,"gibbs-noise")!=0;

    if (!have_latent && !have_locations 
     && model!=0 && (model->type!='R' || model->noise.alpha[2]!=0))
    { fprintf(stderr,
       "Can't sample hyperparameters for a non-Gaussian model when\n");
      fprintf(stderr,
       "both latent vectors and node locations are absent\n");
      exit(1);
    }

    if (!have_locations && dft->N_trees>1)
    { fprintf(stderr,
        "Can't sample to fill in locations when there's more than one tree\n"
      );
      exit(1);
    }

    had_locations = have_locations;
    had_latent = have_latent;
    
    if (!have_locations)
    { sample_locations();
    }
    if (!have_latent) 
    { sample_latent();
    }

    dft_sample_terminals (dft, st, latent, terminals);

    K = (int)pm<=0 ? 1 : (int)pm;

    for (k = 0; k<K; k++)
    { if (do_hypers)
      { gibbs_d_SD();
      }
      if (model!=0 && model->type=='R' && do_noise) 
      { gibbs_noise();
      }
    }

    have_locations = had_locations;
    have_latent = had_latent;

    return 1;
  }

  else if (strcmp(op,"slice-div")==0)
  {
    if (pm<0) 
    { fprintf(stderr,"First parameter of slice-div must be positive\n");
      exit(1);
    }

    slice_div (it, pm==0 ? 1 : pm, (int)pm2<=0 ? 1 : (int)pm2);
    return 1;
  }

  else if (strcmp(op,"sample-latent")==0)
  {
    sample_latent();
    return 1;
  }

  else if (strcmp(op,"sample-locations")==0)
  {
    sample_locations();
    return 1;
  }

  else if (strcmp(op,"gibbs-latent")==0)
  {
    if (have_locations)
    { dft_sample_latent (dft, model, st, inv_temp, latent);
      have_latent = 1;
    }
    else
    { if (dft->N_trees>1)
      { fprintf(stderr,
"Can't do gibbs-latent with no node locations if there's more than one tree\n");
        exit(1);
      }
      if (!have_latent)
      { create_latent();
      }

      K = (int)pm<=0 ? 1 : (int)pm;

      for (k = 0; k<K; k++)
      { dft_gibbs_latent (dft, model, st, inv_temp, latent);
      }
    }

    return 1;
  }

  else if (strcmp(op,"gibbs-locations")==0)
  {
    if (!have_locations)
    { create_locations();
      dft_setup_state (dft, model, st, hyp, parents, divt, 
                       locations, nodes, N_train);
    }

    K = (int)pm<=0 ? 1 : (int)pm;

    for (k = 0; k<K; k++)
    { dft_gibbs_locations (dft, model, st, inv_temp, have_latent ? latent : 0);
    }

    return 1;
  }

  else if (strcmp(op,"mh-latent")==0)
  {
    if (pm<0 || pm>1) 
    { fprintf(stderr,"First parameter of mh-latent must be in (0,1]\n");
      exit(1);
    }

    if (dft->N_trees>1)
    { fprintf(stderr,"Can't do mh-latent if there's more than one tree\n");
      exit(1);
    }

    if (have_locations || model==0) 
    { return 1;
    }

    if (!have_latent)
    { create_latent();
    }

    mh_latent (it, pm==0 ? 1 : pm, (int)pm2<=0 ? 1 : (int)pm2);

    return 1;
  }

  else if (strcmp(op,"discard-latent")==0)
  {
    have_latent = 0;
    return 1;
  }

  else if (strcmp(op,"discard-locations")==0)
  {
    have_locations = 0;
    return 1;
  }

  else if (strcmp(op,"create-latent")==0)
  {  
    create_latent();
    return 1;
  }

  else if (strcmp(op,"create-locations")==0)
  {
    create_locations();
    return 1;
  }

  else if (strcmp(op,"slice-positions")==0)
  {
    slice_positions (it, (int)pm==0 ? 1 : (int)pm, 1);
    return 1;
  }

  else if (strcmp(op,"slice-positions2")==0)
  {
    slice_positions (it, (int)pm==0 ? 1 : (int)pm, 2);
    return 1;
  }

  else if (strcmp(op,"met-terminals")==0)
  {
    met_terminals (it, (int)pm==0 ? 1 : (int)pm, 0);
    return 1;
  }

  else if (strcmp(op,"met-nonterminals")==0)
  {
    met_nonterminals (it, (int)pm==0 ? 1 : (int)pm);
    return 1;
  }

  else if (strcmp(op,"met-terminals-uniform")==0)
  {
    met_terminals (it, (int)pm==0 ? 1 : (int)pm, 1);
    return 1;
  }

  return 0;
}


/* EVALUATE POTENTIAL ENERGY AND ITS GRADIENT.  Currently used only to 
   find weights for Annealed Importance Sampling or tempered transitons, 
   for which the gradient is not needed.  The  energy is just minus the 
   log likelihood given latent values, pretending all parameters (including
   the tree structure) are transformed so that their prior is uniform.  
   This works only if there is a data model, and requires that latent values 
   exist. */

void mc_app_energy
( mc_dynamic_state *ds,	/* Current dynamical state */
  int N_approx,		/* Number of gradient approximations in use */
  int w_approx,		/* Which approximation to use this time */
  double *energy,	/* Place to store energy, null if not required */
  mc_value *gr		/* Place to store gradient, null if not required */
)
{
  double alpha, inv_temp, lat, targ, d;
  int i, t;

  inv_temp = !ds->temp_state ? 1 : ds->temp_state->inv_temp;

  if (gr!=0) abort();

  if (model==0)
  { fprintf (stderr,
       "Can't use annealing/tempering when there's no data model\n");
    exit(1);
  }

  if (!have_latent)
  { fprintf (stderr,
       "Need latent values to do annealing/tempering\n");
    exit(1);
  }

  if (energy==0) return;

  *energy = 0;

  if (inv_temp==0) return;

  create_trees();

  dft_setup_state (dft, model, st, hyp, parents, divt, 
                   have_locations ? locations : 0, nodes, N_train);

  for (t = 0; t<N_targets; t++)
  { for (i = 0; i<N_train; i++)
    { 
      lat = latent[i*N_targets+t];
      targ = train_targets[i*N_targets+t];

      switch (model->type)
      {
        case 'R':
        { alpha = model->noise.alpha[2];
          d = (lat-targ)/st->noise[t];
          if (alpha==0)
          { *energy += log(sqrt(2*M_PI)*st->noise[t]) + d*d/2;
          }
          else
          { *energy += lgamma(alpha/2) - lgamma((alpha+1)/2)
                        + log(sqrt(2*M_PI)*st->noise[t])
                        + ((1+alpha)/2) * log(1+d*d/alpha);
          }
          break;
        }

        case 'B':
        { *energy += targ==1 ? log(1+exp(-lat)) : log(1+exp(lat));
          break;
        }

        default: 
        { abort();
        }
      }
    }
  }

  *energy *= inv_temp;
}


/* SAMPLE FROM DISTRIBUTION AT INVERSE TEMPERATURE OF ZERO.  Returns zero
   if this is not possible. */

int mc_app_zero_gen
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
  int first[Max_trees];
  int N, dt, t;

  create_trees();

  dft_setup_state (dft, model, st, hyp, parents, divt, 
                   locations, nodes, N_train);

  dft_sample_hyper_prior (dft, model, hyp);

  for (N = 0; N<N_train; N++)
  {
    /* Initialize latent vector for new case to zero. */

    for (t = 0; t<N_targets; t++)
    { latent[N*N_targets+t] = 0; 
    }

    /* Add one more case to each of the trees, incrementing latent vector. */

    for (dt = 0; dt<dft->N_trees; dt++)
    { 
      double *vp, *vc, *vn, *v, *w;
      double s;
      int x;

      v = locations + (N_train-1)*N_targets*dt;
      w = terminals[dt];

      if (N==0)
      {
        first[dt] = 1;
        st[dt].parents[1] = 0;

        for (t = 0; t<N_targets; t++)
        { w[N*N_targets+t] = st[dt].d_SD[t]*rand_gaussian();
        }

        dft_conv_tree(st[dt].parents,st[dt].nodes,1);
      }
      else
      {
        dft_gen_path (hyp, st, dt, first[dt], &x, &s, (int *) 0);
  
        first[dt] = dft_add_node (st[dt].parents, st[dt].nodes, -N, N+1, x);
        st[dt].divt[N] = s;

        vp = st[dt].parents[-N]==0 ? 0 : v + (-st[dt].parents[-N]-1)*N_targets;
        vc = x>0 ? w+(x-1)*N_targets : v+(-x-1)*N_targets;
        vn = v + (N-1)*N_targets;

        for (t = 0; t<N_targets; t++)
        { double am, bm, av, bv, sd;
          sd = st[dt].d_SD[t];
          am = vp==0 ? 0 : vp[t];
          av = (sd*sd) * (vp==0 ? s : s - st[dt].divt[-st[dt].parents[-N]]);
          bm = vc[t];
          bv = (sd*sd) * (x>0 ? 1-s : st[dt].divt[-x] - s);
          vn[t] = av+bv<=0 ? (am+bm)/2  /* check for this just in case */
                : (am*bv+bm*av)/(av+bv) + rand_gaussian()*sqrt((av*bv)/(av+bv));
          w[N*N_targets+t] = vn[t] + sqrt(1-s)*sd*rand_gaussian();
        }
      }
  
      for (t = 0; t<N_targets; t++)
      { latent[N*N_targets+t] += w[N*N_targets+t];
        if (isnan(latent[N*N_targets+t])) abort();
      }
    }
  }

  have_latent = 1;

  return 1;
}


/* SET STEPSIZES FOR EACH COORDINATE.  Nothing to do, since there are
   no coordinates. */

void mc_app_stepsizes
( mc_dynamic_state *ds	/* Current dynamical state */
)
{ 
}


/* CREATE TREE STRUCTURES. */

static void create_trees (void)
{
  int dt, i;

  if (N_train!=0 && !have_trees)
  {
    for (i = 0; i<dft->N_trees*(2*N_train); i++)
    { parents[i] = 0;
    }
  
    for (dt = 0; dt<dft->N_trees; dt++)
    { 
      int c, j;
      int *p;
  
      p = parents + (2*N_train)*dt + N_train-1;
        
      for (i = 1; i<=N_train-1; i++)
      { 
        c = rand_int(N_train-i+1);
        for (j = -i+1; ; j++)
        { if (j==0) continue;
          if (j>N_train) abort();
          if (p[j]==0) 
          { c -= 1;
            if (c<0) break;
          }
        }
        p[j] = -i;
  
        c = rand_int(N_train-i);
        for (j = -i+1; ; j++)
        { if (j==0) continue;
          if (j>N_train) abort();
          if (p[j]==0) 
          { c -= 1;
            if (c<0) break;
          }
        }
        p[j] = -i;

        divt [(N_train-1)*dt + i-1] = 0.1 * (double) (N_train-i) / N_train;
      }
    }

    have_trees = 1;
  }
  
  if (have_trees && nodes==0)
  {
    nodes = chk_alloc (dft->N_trees*N_train, sizeof *nodes);
  
    for (dt = 0; dt<dft->N_trees; dt++)
    { (void) dft_conv_tree (parents+dt*(2*N_train)+N_train-1, 
                            nodes+dt*N_train, N_train);
    }
  }
}


/* CREATE LATENT VECTORS. */

static void create_latent (void)
{
  int i, j;
    
  if (have_latent) return;

  for (i = 0; i<dft->N_targets; i++) 
  { 
    for (j = 0; j<N_train; j++)
    { latent[j*dft->N_targets+i] = 
        model!=0 && model->type=='B' ? 2*train_targets[j*dft->N_targets+i] - 1
                                     : train_targets[j*dft->N_targets+i];
    }
  }
  
  have_latent = 1;
}


/* CREATE NODE LOCATIONS. */

static void create_locations (void)
{
  int i, j;

  if (have_locations) return;

  for (i = 0; i<dft->N_targets; i++) 
  { 
    for (j = 0; j<dft->N_trees*(N_train-1); j++)
    { locations[j*dft->N_targets+i] = 0;
    }
  }

  have_locations = 1;
}


/* SAMPLE FOR LATENT VECTORS. */

static void sample_latent (void)
{
  if (!have_locations)
  { 
    if (model!=0 && (model->type!='R' || model->noise.alpha[2]!=0))
    { fprintf(stderr,
"Can't sample for latent vectors with no node locations if model is non-Gaussian\n");
      exit(1);
    }

    if (dft->N_trees>1)
    { fprintf(stderr,
"Can't sample for latent vectors with no node locations with more than one tree\n");
      exit(1);
    }

    dft_setup_state (dft, model, st, hyp, parents, divt, 
                     locations, nodes, N_train);

    dft_gibbs_locations (dft, model, st, inv_temp, 0);
  }
    
  dft_sample_latent (dft, model, st, inv_temp, latent);

  have_latent = 1;
}


/* SAMPLE FOR NODE LOCATIONS. */

static void sample_locations (void)
{
  if (dft->N_trees>1)
  { fprintf(stderr,
      "Can't sample for node locations when there is more than one tree\n");
    exit(1);
  }

  if (!have_locations)
  { have_locations = 1;
    dft_setup_state (dft, model, st, hyp, parents, divt, 
                     locations, nodes, N_train);
  }

  dft_gibbs_locations (dft, model, st, inv_temp, have_latent ? latent : 0);
}


/* DO GIBBS SAMPLING FOR NOISE STANDARD DEVIATIONS. */

static void gibbs_noise (void)
{
  double nalpha, nprec, sum, d, ps;
  prior_spec *pr;
  int i, j;

  pr = &model->noise;

  if (pr->alpha[1]!=0 && pr->alpha[2]==0)
  {
    for (j = 0; j<N_targets; j++)
    {
      sum = pr->alpha[1] * (hyp->noise_cm * hyp->noise_cm);
      for (i = 0; i<N_train; i++)
      { d = latent[i*N_targets+j] - train_targets[i*N_targets+j];
        sum += inv_temp * d*d;
      }

      nalpha = pr->alpha[1] + inv_temp * N_train;
      nprec = nalpha / sum;

      st->noise[j] = prior_pick_sigma (1/sqrt(nprec), nalpha);
    }
  }

  if (pr->alpha[1]!=0 && pr->alpha[2]!=0)
  {
    for (j = 0; j<N_targets; j++)
    {
      ps = pr->alpha[2] * (st->noise[j] * st->noise[j]);

      sum = 0;
      for (i = 0; i<N_train; i++)
      { d = latent[i*N_targets+j] - train_targets[i*N_targets+j];
        sum += rand_gamma((pr->alpha[2]+inv_temp)/2) / ((ps+inv_temp*d*d)/2);
      }

      st->noise[j] = cond_sigma (hyp->noise_cm, pr->alpha[1],
                                    pr->alpha[2], sum, N_train);
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]==0 && pr->alpha[2]==0)
  {
    sum = pr->alpha[0] * (pr->width * pr->width);
    for (i = 0; i<N_train; i++)
    { for (j = 0; j<N_targets; j++)
      { d = latent[i*N_targets+j] - train_targets[i*N_targets+j];
        sum += inv_temp * d*d;
      }
    }

    nalpha = pr->alpha[0] + inv_temp * N_train * N_targets;
    nprec = nalpha / sum;
    hyp->noise_cm = prior_pick_sigma (1/sqrt(nprec), nalpha);

    for (j = 0; j<N_targets; j++)
    { st->noise[j] = hyp->noise_cm;
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]==0 && pr->alpha[2]!=0)
  {
    ps = pr->alpha[2] * (hyp->noise_cm * hyp->noise_cm);

    sum = 0;
    for (i = 0; i<N_train; i++)
    { for (j = 0; j<N_targets; j++) 
      { d = latent[i*N_targets+j] - train_targets[i*N_targets+j];
        sum += rand_gamma((pr->alpha[2]+inv_temp)/2) / ((ps+inv_temp*d*d)/2);
      }
    }

    hyp->noise_cm = cond_sigma (pr->width, pr->alpha[0],
                                   pr->alpha[2], sum, N_targets*N_train);

    for (j = 0; j<N_targets; j++)
    { st->noise[j] = hyp->noise_cm;
    }
  }

  if (pr->alpha[0]!=0 && pr->alpha[1]!=0)
  {
    sum = 0;
    for (j = 0; j<N_targets; j++) 
    { sum += 1 / (st->noise[j] * st->noise[j]);
    }

    hyp->noise_cm = cond_sigma (pr->width, pr->alpha[0],
                                   pr->alpha[1], sum, N_targets);
  }
}


/* DO GIBBS SAMPLING FOR DIFFUSION STANDARD DEVIATIONS. */

static void gibbs_d_SD (void)
{
  double nalpha, nprec, sum, d, s;
  prior_spec *pr;
  int dt, i, x, px;

  for (dt = 0; dt<dft->N_trees; dt++)
  {  
    if (st[dt].locations==0) abort();

    pr = &dft->tree[dt].d_SD;

    if (pr->alpha[1]!=0)
    {
      for (i = 0; i<N_targets; i++)
      {
        sum = pr->alpha[1] * (hyp->d_SD_cm[dt] * hyp->d_SD_cm[dt]);
        nalpha = pr->alpha[1];

        for (x = -(N_train-1); x<=N_train; x++)
        { if (x==0) continue;
          px = st[dt].parents[x];
          if (x>0)
          { d = terminals[dt][N_targets*(x-1)+i] 
                 - (px==0 ? 0 : st[dt].locations[N_targets*(-px-1)+i]);
            s = 1 - (px==0 ? 0 : st[dt].divt[-px]);
          }
          else
          { d = st[dt].locations[N_targets*(-x-1)+i]
                 - (px==0 ? 0 : st[dt].locations[N_targets*(-px-1)+i]);
            s = st[dt].divt[-x] - (px==0 ? 0 : st[dt].divt[-px]);
          }

          if (s>0) /* Disregard singular cases */
          { sum += (d*d) / s;
            nalpha += 1;
          }
        }

        nprec = nalpha / sum;

        st[dt].d_SD[i] = prior_pick_sigma (1/sqrt(nprec), nalpha);
      }
    }
  
    if (pr->alpha[0]!=0 && pr->alpha[1]==0)
    {
      sum = pr->alpha[0] * (pr->width * pr->width);
      nalpha = pr->alpha[0];

      for (i = 0; i<N_targets; i++)
      { for (x = -(N_train-1); x<=N_train; x++)
        { if (x==0) continue;
          px = st[dt].parents[x];
          if (x>0)
          { d = terminals[dt][N_targets*(x-1)+i] 
                 - (px==0 ? 0 : st[dt].locations[N_targets*(-px-1)+i]);
            s = 1 - (px==0 ? 0 : st[dt].divt[-px]);
          }
          else
          { d = st[dt].locations[N_targets*(-x-1)+i]
                 - (px==0 ? 0 : st[dt].locations[N_targets*(-px-1)+i]);
            s = st[dt].divt[-x] - (px==0 ? 0 : st[dt].divt[-px]);
          }
          if (s>0) /* Disregard singular cases */
          { sum += (d*d) / s;
            nalpha += 1;
          }
        }
      }

      nprec = nalpha / sum;
  
      hyp->d_SD_cm[dt] = prior_pick_sigma (1/sqrt(nprec), nalpha);
  
      for (i = 0; i<N_targets; i++)
      { st[dt].d_SD[i] = hyp->d_SD_cm[dt];
      }
    }
  
    if (pr->alpha[0]!=0 && pr->alpha[1]!=0)
    {
      sum = 0;
      for (i = 0; i<N_targets; i++) 
      { sum += 1 / (st[dt].d_SD[i] * st[dt].d_SD[i]);
      }

      hyp->d_SD_cm[dt] 
        = cond_sigma (pr->width, pr->alpha[0], pr->alpha[1], sum, N_targets);
    }
  }
}


/* DO SLICE SAMPLING FOR PARAMETERS OF DIVERGENCE FUNCTION. */

static void slice_div
( mc_iter *it,		/* Description of this iteration */
  double w,		/* Width of initial interval */
  int K			/* Number of scans to do */
)
{
  double cur_ll, new_ll, slice_lp, new_lp, cur_val, new_val, low_val, high_val;
  prior_spec *prior;
  double *param, alpha, omega;
  int dt, k, c;
  int first;

  for (dt = 0; dt<dft->N_trees; dt++)
  { 
    first = 1;

    for (k = 0; k<K; k++)
    { for (c = 0; c<3; c++)
      { 
        prior = c==0 ? &dft->tree[dt].c0 
              : c==1 ? &dft->tree[dt].c1
              :        &dft->tree[dt].c2;

        alpha = prior->alpha[0];

        if (prior->width==0 || alpha==0) 
        { continue;
        }

        omega = 1/(prior->width*prior->width);

        if (first)
        { 
          cur_ll = dft_log_prob_div (dft, hyp, dt, st, 0, 0);

          /* Refuse to continue with sampling if current log prob for
             divergences is -INFINITY, due to divergence at time 1.  
             It shouldn't ever be +INFINITY. */

          if (cur_ll==-INFINITY) return;
          if (cur_ll==INFINITY) abort();

          it->slice_evals += 1;
          first = 0;
        }

        param = c==0 ? &hyp->c0[dt] 
              : c==1 ? &hyp->c1[dt] 
              :        &hyp->c2[dt];

        it->slice_calls += 1;

        cur_val = -2*log(*param);
        low_val = cur_val - w*rand_uniopen();
        high_val = low_val + w;

        slice_lp = cur_val*alpha/2 - exp(cur_val)*alpha/(2*omega)
                 + cur_ll - rand_exp();

        for (;;)
        { 
          new_val = low_val + (high_val-low_val)*rand_uniopen();
          *param = exp(-new_val/2);

          new_ll = dft_log_prob_div (dft, hyp, dt, st, 0, 0);
          new_lp = new_val*alpha/2 - exp(new_val)*alpha/(2*omega) + new_ll;
          it->slice_evals += 1;

          if (new_lp>=slice_lp) 
          { break;
          }

          if (new_val>cur_val)
          { high_val = new_val;
          }
          else
          { low_val = new_val;
          }
        }

        cur_ll = new_ll;
      }
    }
  }
}


/* METROPOLIS-HASTINGS UPDATES FOR LATENT VECTORS. */

static double latent_likelihood
( double target,
  double value,
  int t
)
{ double lk;
  
  switch (model->type)
  { case 'R':
    { double d, a;
      d = (target-value)*(target-value) / (st->noise[t]*st->noise[t]);
      a = model->noise.alpha[2];
      lk = a==0 ? -0.5*d : -0.5*(a+1)*log(1+d/a);
      break;
    }
    case 'B':
    { lk = target==1 ? -log(1+exp(-value)) : -log(1+exp(value));
      break;
    }
    default:
    { abort();
    }
  }

  return inv_temp * lk;
}

static void mh_latent
( mc_iter *it,		/* Description of this iteration */
  double scale,		/* Scale factor */
  int K			/* Number of scans to do */
)
{ 
  int N_targets;
  double *tmp;
  int i, t, k;
  double f, v, lk0, lk1;
  double *p, *q, *r;

  N_targets = dft->N_targets;
  f = sqrt(1-scale*scale);

  tmp = chk_alloc (N_train*N_targets, sizeof *tmp);

  dft_setup_state (dft, model, st, hyp, parents, divt, 
                   locations, nodes, N_train);

  lk0 = lk1 = 0;

  for (k = 0; k<K; k++)
  {
    dft_sample_prior (dft, st, tmp);

    p = tmp;
    q = latent;
    r = train_targets;

    for (i = 0; i<N_train; i++)
    { for (t = 0; t<N_targets; t++)
      { v = *p;
        *p = *q;
        *q = *q * f + v * scale;
        lk0 += latent_likelihood(*r,*p,t);
        lk1 += latent_likelihood(*r,*q,t);
        p += 1;
        q += 1;
        r += 1;
      }
    }

    it->proposals += 1;
    it->delta = lk0 - lk1;

    if (rand_uniform()<exp(-it->delta))
    { it->move_point = 1;
    }
    else
    { it->move_point = 0;
      it->rejects += 1;
      memcpy (latent, tmp, N_train * N_targets * sizeof(double));
    }
  }

  free(tmp);
}


/* SLICE SAMPLING UPDATE FOR NODE POSITIONS. */

static void slice_positions
( mc_iter *it,		/* Description of this iteration */
  int K,		/* Number of scans to do */
  int method		/* Method for selecting node (1 or 2) */
)
{
  double cur, new, low, high, cur_lev, new_lev, slice_lev, diff;
  int cur_infp, slice_infp, new_infp;
  int dt, x, y, a, b, c, d, k, h;
  dft_likelihood *lka;

  /* Check that required parts are present. */

  if (!have_locations && dft->N_trees>1)
  { fprintf(stderr,
"Can't do slice-positions without node locations if there's more than one tree\n");
    exit(1);
  }
  if (!have_latent 
   && model!=0 && (model->type!='R' || model->noise.alpha[2]!=0))
  { fprintf(stderr,
"Can't do slice-positions without latent vectors if the model is non-Gaussian\n");
    exit(1);
  }

  /* Must have at least two training cases for this to make sense. */

  if (N_train<2) return;

  /* Set up flag array if necessary. */

  if (flags_n!=N_train+1)
  { if (flags!=0) free(flags);
    flags = chk_alloc (N_train+1, sizeof *flags);
    flags_n = N_train+1;
  }

  /* Allocate likelihood array if necessary. */

  if (K>0 && !have_locations)
  { lka = chk_alloc (2*N_train*dft->N_targets, sizeof *lka);
  }
  else
  { lka = 0;
  }

  /* Do K scans, updating each tree in turn. */

  if (K<0) K = -K;

  for (k = 0; k<K; k++)
  { for (dt = 0; dt<dft->N_trees; dt++)
    {
      /* Compute likelihoods non-incrementally whenever we start a new tree. */

      if (lka!=0 && (k==0 && dt==0 || dft->N_trees>1))
      { double mean, var, peak;
        int t, infp;
        for (t = 0; t<dft->N_targets; t++)
        { dft_tree_likelihood (dft, dt, st, have_latent ? latent : 0,
                               t, 0, inv_temp, &mean, &var, &peak, &infp,
                               lka + 2*N_train*t + N_train-1);
        }
      }

      /* Do updates w.r.t. paths to each terminal node, y, in turn. */
  
      for (y = 1; y<=N_train; y++)
      { 
        it->slice_calls += 1;

        /* Find node along path to update. */

        switch (method)
        { 
          case 1: /* Random position along path */
          {
            int cnt;

            cnt = 0;
            for (a = st[dt].parents[y]; a!=0; a = st[dt].parents[a]) 
            { cnt += 1;
            }

            b = y;
            a = st[dt].parents[b]; 
            for (cnt = rand_int(cnt); cnt>0; cnt--)
            { if (a==0) abort();
              b = a;
              a = st[dt].parents[b];
            }

            b = dft_sibling (st[dt].parents, st[dt].nodes, b);

            break;
          }

          case 2: /* Common anscestor with other terminal node */
          {
            /* Pick another terminal node, x. */

            x = 1 + rand_int(N_train-1);
            if (x>=y) x += 1;
    
            /* Find the most recent ancestor, a, of x and y, along with the
               child, b, of a on the path to x (note that b may be x). */
    
            a = st[dt].parents[x];
            while (a!=0)
            { flags[-a] = 0;
              a = st[dt].parents[a];
            }
            a = st[dt].parents[y];
            while (a!=0)
            { flags[-a] = 1;
              a = st[dt].parents[a];
            }
            a = st[dt].parents[x];
            b = x;
            while (flags[-a]==0)
            { b = a;
              a = st[dt].parents[a];
              if (a==0) abort();
            }

            break;
          }

          default:
          { abort();
          }
        }

        /* Find the slice level based on current position. */

        cur = st[dt].divt[-a];
        if (cur==1.0 && dft_div(hyp,dt,1.0)>0) 
        { continue; /* Point mass at t=1, which isn't in (0,1). */
        }

        dft_log_prob_node (dft, dt, st, have_latent ? latent : 0, b, 
                           inv_temp, &cur_lev, &cur_infp, lka);
        if (cur_lev==-INFINITY) abort();
        cur_lev += dft_log_prob_paths (dft, hyp, dt, st, b);
        if (cur_lev==-INFINITY) abort();

        it->slice_evals += 1;

        slice_lev = cur_lev - rand_exp();
        slice_infp = cur_infp;

        /* Remove the subtree off a headed by b. */
   
        if (lka!=0)
        { c = dft_sibling(st[dt].parents,st[dt].nodes,b);
        }
        dft_remove_node (st[dt].parents, st[dt].nodes, a, b);
        if (lka!=0)
        { dft_update_likelihoods (dft, dt, st, have_latent ? latent : 0, 
                                  st[dt].parents[c], lka);
        }

        /* Set up initial interval. */

        low = 0;
        high = b>0 ? 1 : st[dt].divt[-b];
        h = y;

        /* Repeatedly sample from the interval, and either accept the 
           point drawn, or use it to narrow the interval. */

        for (;;)
        { 
          /* Draw a new divergence time from current interval. */

          new = low + (high-low)*rand_uniopen();

          if (new<=low || new>=high) /* ensure against problems from roundoff */
          { new = (low+high)/2;
          }

          /* Find the node within the path to y that is at the end of the
             segment containing this divergence time. */

          c = h;
          for (;;)
          { d = st[dt].parents[c];
            if (d==0 || st[dt].divt[-d]<=new) break;
            c = d;
          }

          /* Add a with subtree headed by b back into tree at this position. */

          dft_add_node (st[dt].parents, st[dt].nodes, a, b, c);
          st[dt].divt[-a] = new;
          if (lka!=0)
          { dft_update_likelihoods (dft, dt, st, have_latent ? latent : 0, 
                                    a, lka);
          }

          /* Find probability level with a in new position. */

          dft_log_prob_node (dft, dt, st, have_latent ? latent : 0, b, 
                             inv_temp, &new_lev, &new_infp, lka);
          new_lev += dft_log_prob_paths (dft, hyp, dt, st, b);

          it->slice_evals += 1;

          /* Accept if level is at least as great as slice level. */

          if (new_lev!=-INFINITY && (new_infp>slice_infp 
                || new_infp==slice_infp && new_lev>=slice_lev))
          { break;
          }

          /* If point is rejected, remove a again, and narrow interval. */

          dft_remove_node (st[dt].parents, st[dt].nodes, a, b);
          if (lka!=0 && st[dt].parents[c]!=0)
          { dft_update_likelihoods (dft, dt, st, have_latent ? latent : 0, 
                                    st[dt].parents[c], lka);
          }

          if (new>cur)
          { high = new;
            h = c;
          }
          else
          { low = new;
          }
        }
      }
    }
  }

  if (lka) free(lka);
}


/* METROPOLIS UPDATES FOR TERMINAL NODES.  If uniform is zero, the new 
   position of the parent of the node being updated is chosen by simulating 
   the generation process; if uniform is one, the new position is choosen 
   by picking a segment uniformly and then a time uniformly. */

static void met_terminals
( mc_iter *it,		/* Description of this iteration */
  int K,		/* Number of scans to do */
  int uniform		/* Choose attachment point uniformly? */
)
{
  double logprob0, logprob1;
  double odivt, t1, t0;
  int infp0, infp1;
  int dt, x, y, yc0, yc1, k;
  dft_likelihood *lka;
  int first;

  /* Check that required parts are present. */

  if (!have_locations && dft->N_trees>1)
  { fprintf(stderr,
"Can't do met-terminals without node locations if there's more than one tree\n");
    exit(1);
  }
  if (!have_latent 
   && model!=0 && (model->type!='R' || model->noise.alpha[2]!=0))
  { fprintf(stderr,
"Can't do met-terminals without latent vectors if the model is non-Gaussian\n");
    exit(1);
  }

  /* Must have at least two training cases for this to make sense. */

  if (N_train<2) return;

  /* Allocate likelihood array if necessary. */

  if (K>0 && !have_locations)
  { lka = chk_alloc (2*N_train*dft->N_targets, sizeof *lka);
  }
  else
  { lka = 0;
  }

  /* Do K scans, updating each tree in turn. */

  if (K<0) K = -K;

  for (k = 0; k<K; k++)
  { for (dt = 0; dt<dft->N_trees; dt++)
    {
      /* Compute likelihoods non-incrementally whenever we start a new tree. */

      if (lka!=0 && (k==0 && dt==0 || dft->N_trees>1))
      { double mean, var, peak;
        int t, infp;
        for (t = 0; t<dft->N_targets; t++)
        { dft_tree_likelihood (dft, dt, st, have_latent ? latent : 0,
                               t, 0, inv_temp, &mean, &var, &peak, &infp,
                               lka + 2*N_train*t + N_train-1);
        }
      }

      /* Do updates for each terminal node, x, in turn. */
  
      for (x = 1; x<=N_train; x++)
      { 
        /* Find the parent, y, of x, and its current other child, yc0. */
  
        y = st[dt].parents[x];
        yc0 = dft_sibling (st[dt].parents, st[dt].nodes, x);

        if (uniform && st[dt].divt[-y]==1.0 && dft_div(hyp,dt,1.0)>0) 
        { continue; /* Point mass at t=1, which lacks +ve proposal probability*/
        }
  
        /* Compute log probability of generation of x, omitting the probability 
           for the generation of the path in the tree to x if we are
           choosing the new position by simulating this.  Also generates a 
           new location for the parent of x, if locations are present. */
  
        dft_log_prob_node (dft, dt, st, have_latent ? latent : 0, x, 
                           inv_temp, &logprob0, &infp0, lka);
        if (logprob0==-INFINITY) abort();
        if (uniform) 
        { logprob0 += dft_log_prob_path (dft, hyp, dt, st, x);
          if (logprob0==-INFINITY) abort();
        }
  
        /* If we are choosing the new position uniformly, adjust logprob0 to 
           account for asymmetrical divergence time proposals, by adding the 
           log of the time interval for the segment. */
  
        if (uniform)
        { t0 = st[dt].parents[y]==0 ? 0 : st[dt].divt[-st[dt].parents[y]];
          t1 = yc0>0 ? 1 : st[dt].divt[-yc0];
          if (t0>t1 || t1>1 || t0<0) abort();
          logprob0 += log(t1-t0);
        }
  
        /* Save current divergence time and location of parent of x. */
  
        odivt = st[dt].divt[-y];
  
        if (st[dt].locations)
        { memcpy (saved_location, st[dt].locations + (-y-1)*dft->N_targets, 
                  dft->N_targets * sizeof (double));
        }
  
        /* Remove y (with x) from the tree. */
  
        first = dft_remove_node (st[dt].parents, st[dt].nodes, y, x);
        if (lka!=0)
        { dft_update_likelihoods (dft, dt, st, have_latent ? latent : 0, 
                                  st[dt].parents[yc0], lka);
        }
  
        /* Randomly pick the node, yc1, to become the new other child of y. 
           Could turn out to be the same as the old other child.  This
           either done uniformly or by simulating the generation process. 
           Also randomly pick the new divergence time, and adjust logprob0
           if necessary to account for asymmetric proposals. */
  
        if (uniform)
        { 
          double new;

          do 
          { yc1 = rand_int(2*N_train-1) - (N_train-1);
            if (yc1>=0) yc1 += 1;
          } while (yc1==x || yc1==y);
  
          t0 = st[dt].parents[yc1]==0 ? 0 : st[dt].divt[-st[dt].parents[yc1]];
          t1 = yc1>0 ? 1 : st[dt].divt[-yc1];
          if (t0>t1 || t1>1 || t0<0) abort();
  
          new = t0 + (t1-t0)*rand_uniopen();
          if (new<=t0 || new>=t1) /* ensure against problems from roundoff */
          { new = (t0+t1)/2;
          }
          st[dt].divt[-y] = new;
  
          logprob0 -= log(t1-t0);
        }
  
        else /* not uniform */
        { 
          dft_gen_path (hyp, st, dt, first, &yc1, &st[dt].divt[-y], (int *) 0);
        }

        /* Add y (with x) back to the tree as parent of yc1. */
  
        dft_add_node (st[dt].parents, st[dt].nodes, y, x, yc1);
        if (lka!=0)
        { dft_update_likelihoods (dft, dt, st, have_latent ? latent : 0, 
                                  y, lka);
        }
  
        /* Compute new log probability of generation of x.  Also randomly 
           generates a new location for the parent of x, if locations are 
           present. */
  
        dft_log_prob_node (dft, dt, st, have_latent ? latent : 0, x, 
                           inv_temp, &logprob1, &infp1, lka);
        if (uniform) 
        { logprob1 += dft_log_prob_path (dft, hyp, dt, st, x);
        }

        /* Decide whether to reject, and if so restore old tree. */
  
        it->proposals += 1;
        it->delta = (logprob0 - logprob1) + 1e10 * (infp0 - infp1);

        if (logprob1!=-INFINITY 
         && (infp1>infp0 || infp1==infp0 && rand_uniform()<exp(-it->delta)))
        { it->move_point = 1;
        }
        else
        { dft_remove_node (st[dt].parents, st[dt].nodes, y, x);
          if (lka!=0)
          { dft_update_likelihoods(dft, dt, st, have_latent ? latent : 0, 
                                   st[dt].parents[yc1], lka);
          }
          dft_add_node (st[dt].parents, st[dt].nodes, y, x, yc0);
          st[dt].divt[-y] = odivt;
          if (lka!=0)
          { dft_update_likelihoods(dft, dt, st, have_latent ? latent : 0, 
                                   y, lka);
          }
          if (st[dt].locations)
          { memcpy (st[dt].locations + (-y-1)*dft->N_targets, saved_location,
                    dft->N_targets * sizeof (double));
          }
          it->move_point = 0;
          it->rejects += 1;
        }
      }
    }
  }

  if (lka) free(lka);
}


/* METROPOLIS UPDATES FOR NONTERMINAL NODES. */

static void met_nonterminals
( mc_iter *it,		/* Description of this iteration */
  int K			/* Number of scans to do */
)
{
  double logprob0, logprob1;
  double odivt, t1, t0;
  int infp0, infp1;
  int dt, x, y, yc0, yc1, k;
  dft_likelihood *lka;
  int first, cnt;

  /* Check that required parts are present. */

  if (!have_locations && dft->N_trees>1)
  { fprintf(stderr,
"Can't do met-nonterminals without node locations if there's more than one tree\n");
    exit(1);
  }
  if (!have_latent 
   && model!=0 && (model->type!='R' || model->noise.alpha[2]!=0))
  { fprintf(stderr,
"Can't do met-nonterminals without latent vectors if the model is non-Gaussian\n");
    exit(1);
  }

  /* Must have at least two training cases for this to make sense. */

  if (N_train<2) return;

  /* Allocate likelihood array if necessary. */

  if (K>0 && !have_locations)
  { lka = chk_alloc (2*N_train*dft->N_targets, sizeof *lka);
  }
  else
  { lka = 0;
  }

  /* Do K scans, updating each tree in turn. */

  if (K<0) K = -K;

  for (k = 0; k<K; k++)
  { for (dt = 0; dt<dft->N_trees; dt++)
    {
      /* Compute likelihoods non-incrementally whenever we start a new tree. */

      if (lka!=0 && (k==0 && dt==0 || dft->N_trees>1))
      { double mean, var, peak;
        int t, infp;
        for (t = 0; t<dft->N_targets; t++)
        { dft_tree_likelihood (dft, dt, st, have_latent ? latent : 0,
                               t, 0, inv_temp, &mean, &var, &peak, &infp,
                               lka + 2*N_train*t + N_train-1);
        }
      }

      /* Do updates for each nonterminal node, x, in turn, except the root. */
  
      for (x = -1; x>=-(N_train-1); x--)
      { 
        if (st[dt].parents[x]==0) continue;

        /* Find the parent, y, of x, and its current other child, yc0. */
  
        y = st[dt].parents[x];
        yc0 = dft_sibling (st[dt].parents, st[dt].nodes, x);

        /* Compute log probability of generation of the subtree headed by x, 
           omitting the probability for the generation of the first path in 
           the tree to a node in x, and the probabilities internal to the
           subtree headed by x.  Also generates a new location for the parent 
           of x, if locations are present. */
  
        dft_log_prob_node (dft, dt, st, have_latent ? latent : 0, x, 
                           inv_temp, &logprob0, &infp0, lka);
        if (logprob0==-INFINITY) abort();
        logprob0 += dft_log_prob_paths (dft, hyp, dt, st, x);
        logprob0 -= dft_log_prob_path (dft, hyp, dt, st, x);
        if (logprob0==-INFINITY) abort();
  
        /* Save current divergence time and location of parent of x. */
  
        odivt = st[dt].divt[-y];
  
        if (st[dt].locations)
        { memcpy (saved_location, st[dt].locations + (-y-1)*dft->N_targets, 
                  dft->N_targets * sizeof (double));
        }
  
        /* Remove y (with x) from the tree.  Wait till later to update
           likelihoods, since we may just put it back in. */
  
        first = dft_remove_node (st[dt].parents, st[dt].nodes, y, x);
  
        /* Randomly pick the node, yc1, to become the new other child of y. 
           Could turn out to be the same as the old other child.  This is done
           by simulating the generation process.  Also randomly pick the new 
           divergence time.  This process is repeated if the new divergence
           time doesn't precede the divergence time for x, up to a maximum
           of 100 times, which if reached produces rejection. */

        for (cnt = 0; cnt<100; cnt++)
        { 
          dft_gen_path (hyp, st, dt, first, &yc1, &st[dt].divt[-y], (int *) 0);
          if (st[dt].divt[-y]<st[dt].divt[-x]) break;
        }

        if (cnt==100) 
        { 
          dft_add_node (st[dt].parents, st[dt].nodes, y, x, yc0);
          st[dt].divt[-y] = odivt;
          it->move_point = 0;
          it->proposals += 1;
          it->rejects += 1;

          continue; /* go on to next x */
        }

        /* Add y (with x) back to the tree as parent of yc1, and update
           likelihoods. */
  
        if (lka!=0)
        { dft_update_likelihoods (dft, dt, st, have_latent ? latent : 0, 
                                  st[dt].parents[yc0], lka);
        }
        dft_add_node (st[dt].parents, st[dt].nodes, y, x, yc1);
        if (lka!=0)
        { dft_update_likelihoods (dft, dt, st, have_latent ? latent : 0, 
                                  y, lka);
        }
  
        /* Compute new log probability of generation of x.  Also randomly 
           generates a new location for the parent of x, if locations are 
           present. */
  
        dft_log_prob_node (dft, dt, st, have_latent ? latent : 0, x, 
                           inv_temp, &logprob1, &infp1, lka);
        logprob1 += dft_log_prob_paths (dft, hyp, dt, st, x);
        logprob1 -= dft_log_prob_path (dft, hyp, dt, st, x);

        /* Decide whether to reject, and if so restore old tree. */
  
        it->proposals += 1;
        it->delta = (logprob0 - logprob1) + 1e10 * (infp0 - infp1);

        if (logprob1!=-INFINITY 
         && (infp1>infp0 || infp1==infp0 && rand_uniform()<exp(-it->delta)))
        { it->move_point = 1;
        }
        else
        { dft_remove_node (st[dt].parents, st[dt].nodes, y, x);
          if (lka!=0)
          { dft_update_likelihoods(dft, dt, st, have_latent ? latent : 0, 
                                   st[dt].parents[yc1], lka);
          }
          dft_add_node (st[dt].parents, st[dt].nodes, y, x, yc0);
          st[dt].divt[-y] = odivt;
          if (lka!=0)
          { dft_update_likelihoods(dft, dt, st, have_latent ? latent : 0, 
                                   y, lka);
          }
          if (st[dt].locations)
          { memcpy (st[dt].locations + (-y-1)*dft->N_targets, saved_location,
                    dft->N_targets * sizeof (double));
          }
          it->move_point = 0;
          it->rejects += 1;
        }
      }
    }
  }

  if (lka) free(lka);
}
