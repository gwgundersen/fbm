/* DFT-UTIL.C - Utility routines for diffusion tree programs. */

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
#include "data.h"
#include "prior.h"
#include "model.h"
#include "dft.h"


/* SET UP REQUIRED RECORD SIZES.  Doesn't set the sizes of the variable-sized
   records. */

void dft_record_sizes
( log_gobbled *logg     /* Structure to hold gobbled data */
)
{
  logg->req_size['P'] = sizeof (dft_spec);
  logg->req_size['M'] = sizeof (model_specification);
}


/* REPORT ERROR IF DIFFUSION TREE MODEL SPECS, DATA MODEL, ETC. ARE MISSING. */

void dft_check_specs_present
( dft_spec *dft,		/* Diffusion tree model specs, or null */
  int need_model,		/* Must data model be present? */
  model_specification *m	/* Data model, or null */
)
{
  if (dft==0)
  { fprintf(stderr,"No specification for diffusion tree model in log file\n");
    exit(1);
  }

  if (m==0 && need_model)
  { fprintf(stderr,"No model specification in log file\n");
    exit(1);
  }
}


/* COMPUTE SIZE OF HYPERPARAMETERS/PARAMETERS RECORD. */

int dft_hyper_size
( dft_spec *dft,		/* Diffusion tree model specification */
  model_specification *m	/* Data model, or null */
)
{
  int hs;

  hs = sizeof (dft_hypers) + dft->N_targets*dft->N_trees * sizeof (double);

  if (m!=0 && m->type=='R')
  { hs += dft->N_targets * sizeof (double);
  }

  hs -= sizeof (double);  /* Because SD already has one double in the struct */

  return hs;
}


/* COMPUTE SIZE OF DIVERGENCE TIME RECORD. */

int dft_divt_size
( dft_spec *dft,		/* Diffusion tree model specification */
  int N_cases			/* Number of cases */
)
{ return N_cases==0 ? 0 : dft->N_trees * (N_cases-1) * sizeof(double);
}


/* COMPUTE SIZE OF RECORD CONTAINING PARENTS OF NODES.  Has surplus
   entries for non-existent node 0, in order to simplify indexing. */

int dft_parents_size
( dft_spec *dft,		/* Diffusion tree model specification */
  int N_cases			/* Number of cases */
)
{ return dft->N_trees * (2*N_cases) * sizeof(int);
}


/* COMPUTE SIZE OF LATENT VECTOR RECORD. */

int dft_latent_size
( dft_spec *dft,		/* Diffusion tree model specification */
  int N_cases			/* Number of cases */
)
{ return N_cases * dft->N_targets * sizeof(double);
}


/* COMPUTE SIZE OF NODE LOCATIONS RECORD. */

int dft_locations_size
( dft_spec *dft,		/* Diffusion tree model specification */
  int N_cases			/* Number of cases */
)
{ return N_cases==0 ? 0  
          : dft->N_trees * (N_cases-1) * dft->N_targets * sizeof(double);
}


/* CHECK THAT RECORDS HAVE THE RIGHT SIZES. */

void dft_check_sizes
( log_gobbled *logg,		/* Records gobbled from log file */
  dft_spec *dft,		/* Diffusion tree model specification */
  model_specification *model,	/* Data model, or null */
  int N_train			/* Number of training cases */
)
{
  if (logg->data['S']!=0)
  { if (logg->actual_size['S'] != dft_hyper_size(dft,model))
    { fprintf(stderr,
              "Diffusion tree hyperparameter record is the wrong size!\n");
      exit(1);
    }
  }

  if (logg->data['T']!=0)
  { if (logg->actual_size['T'] != dft_divt_size(dft,N_train))
    { fprintf(stderr,"Record of divergence times is the wrong size!\n");
      exit(1);
    }
  }
 
  if (logg->data['R']!=0)
  { if (logg->actual_size['R'] != dft_parents_size(dft,N_train))
    { fprintf(stderr,"Record of node parents is the wrong size!\n");
      exit(1);
    }
  }

  if (logg->data['L']!=0)
  { if (logg->actual_size['L'] != dft_latent_size(dft,N_train))
    { fprintf(stderr,"Record of latent vectors is the wrong size!\n");
      exit(1);
    }
  }

  if (logg->data['N'])
  { if (logg->actual_size['N'] != dft_locations_size(dft,N_train))
    { fprintf(stderr,"Record of node locations is the wrong size!\n");
      exit(1);
    }
  }
}


/* SET UP INITIAL VALUES FOR HYPERPARAMETERS.  Sets the hyperparameters
 and parameters to the location parameters of their priors. */

void dft_hyper_init
( dft_spec *dft,		/* Specification for diffusion tree model */
  model_specification *m,	/* Specification for data model */
  dft_hypers *h			/* Hyperparameter structure to initialize */
)
{
  int dt, t, j;

  j = 0;

  for (dt = 0; dt<dft->N_trees; dt++)
  {
    h->d_SD_cm[dt] = dft->tree[dt].d_SD.width;

    for (t = 0; t<dft->N_targets; t++) 
    { h->SD[j++] = h->d_SD_cm[dt];
    }

    h->c0[dt] = dft->tree[dt].c0.width;
    h->c1[dt] = dft->tree[dt].c1.width;
    h->c2[dt] = dft->tree[dt].c2.width;
  }

  if (m!=0 && m->type=='R')
  {
    h->noise_cm = m->noise.width;

    for (t = 0; t<dft->N_targets; t++)
    { h->SD[j++] = h->noise_cm;
    }
  }
}


/* CONVERT TREE TO INTERNAL FORM.  The value returned is the earliest node
   in the tree (zero if there are no cases in the tree). */

int dft_conv_tree 
( int *parents,			/* Indexes of parents, offset for indexing */
  dft_tree_node *nodes,		/* Array of nodes to store new representation */
  int N_cases			/* Number of cases in tree */
)
{ 
  int i, j, k, c, e;

  if (N_cases==0) return 0;

  for (i = 1; i<=N_cases-1; i++)
  { nodes[i].child[0] = nodes[i].child[1] = 0;
    nodes[i].n_points[0] = nodes[i].n_points[1] = 0;    
  }

  e = 0;

  for (i = 1; i<=N_cases; i++)
  { k = i;
    j = parents[k];
    while (j!=0)
    { c = nodes[-j].child[0]==0 || nodes[-j].child[0]==k ? 0 : 1;
      if (nodes[-j].child[c]==0) 
      { nodes[-j].child[c] = k;
      }
      if (nodes[-j].child[c]!=k) abort();
      nodes[-j].n_points[c] += 1;
      k = j;
      j = parents[k];
    }
    if (e!=0 && e!=k) abort();
    e = k;
  }

  for (i = 1; i<=N_cases-1; i++)
  { if (nodes[i].child[0]==0 || nodes[i].child[1]==0
     || nodes[i].n_points[0]== 0 || nodes[i].n_points[1]==0) 
    { abort();
    }
  }

  if (N_cases>1 && nodes[-e].n_points[0]+nodes[-e].n_points[1]!=N_cases) 
  { abort();
  }

  return e;
}


/* ADD A NODE TO A TREE.  One of the new node's children is already 
   in the tree, the other is not.  The new node's parent is the parent
   of the child presently in the tree.  The "parents" array is passed
   as a pointer to the middle, where non-existent node 0 lies.  The
   value returned is the earliest node in the modified tree. */

int dft_add_node
( int *parents,			/* Pointer to parents, offset for indexing */
  dft_tree_node *nodes,		/* Nodes with child pointers */
  int new_node,			/* Index of new node (negative) */
  int new_child,		/* Index of child not presently in the tree */
  int old_child			/* Index of child already in the tree */
)
{ 
  int parent, e, c;

  parent = parents[old_child];
  parents[new_node] = parent;

  nodes[-new_node].child[0] = new_child;
  nodes[-new_node].child[1] = old_child;

  parents[new_child] = new_node;
  parents[old_child] = new_node;

  if (parent!=0)
  {
    if (nodes[-parent].child[0]==old_child)
    { nodes[-parent].child[0] = new_node;
    }
    else if (nodes[-parent].child[1]==old_child)
    { nodes[-parent].child[1] = new_node;
    }
    else
    { abort();
    }
  }

  e = new_node;
  for (;;)
  { for (c = 0; c<=1; c++)
    { nodes[-e].n_points[c] = 
        nodes[-e].child[c]>0 ? 1 : nodes[-nodes[-e].child[c]].n_points[0] 
                                    + nodes[-nodes[-e].child[c]].n_points[1];
    }
    if (parents[e]==0) break;
    e = parents[e];
  }

  return e;
}


/* REMOVE NODE FROM TREE.  The nonterminal node rm_node is removed from the
   tree, along with its child rm_child.  The other child of rm_node stays
   in the tree, with its new parent being its current grandparent.  The
   value returned is the earliest node in the modified tree. */

int dft_remove_node
( int *parents,			/* Pointer to parents, offset for indexing */
  dft_tree_node *nodes,		/* Nodes with child pointers */
  int rm_node,			/* Index of node to remove (negative) */
  int rm_child			/* Index of child to remove with rm_node */
)
{ 
  int parent, other_child;
  int e, c;

  other_child = dft_sibling (parents, nodes, rm_child);

  parent = parents[rm_node];
  parents[other_child] = parent;

  if (parent!=0)
  {
    if (nodes[-parent].child[0]==rm_node)
    { nodes[-parent].child[0] = other_child;
    }
    else if (nodes[-parent].child[1]==rm_node)
    { nodes[-parent].child[1] = other_child;
    }
    else
    { abort();
    }
  }
  
  e = other_child;

  while (parents[e]!=0)
  { e = parents[e];
    for (c = 0; c<=1; c++)
    { nodes[-e].n_points[c] = 
        nodes[-e].child[c]>0 ? 1 : nodes[-nodes[-e].child[c]].n_points[0] 
                                    + nodes[-nodes[-e].child[c]].n_points[1];
    }
  }

  return e;
}


/* FIND SIBLING OF A NODE. */

int dft_sibling
( int *parents,			/* Pointer to parents, offset for indexing */
  dft_tree_node *nodes,		/* Nodes with child pointers */
  int x				/* Index of node to find sibling of (negative)*/
)
{
  int p;

  if (x==0) abort();

  p = parents[x];

  if (p==0) abort();

  if (nodes[-p].child[0]==x)
  { return nodes[-p].child[1];
  }
  else if (nodes[-p].child[1]==x)
  { return nodes[-p].child[0];
  }
  else
  { abort();
  }
}


/* FIND THE MAXIMUM DEPTH OF A TREE. */

int dft_max_depth
( int *parents,
  int N_cases
)
{ 
  int i, j, d, c;

  d = 0;

  for (i = 1; i<=N_cases; i++)
  { c = 0;
    for (j = i; j!=0; j = parents[j])
    { c += 1;
    }
    if (c>d) d = c;
  }

  return d;
}


/* SET UP POINTERS TO STATE.  Any part may be null, in which case the pointers
   to it are set to be null. */

void dft_setup_state
( dft_spec *dft,		/* Specification of diffusion tree model */
  model_specification *model,	/* Specification of data model */
  dft_state st,			/* Array of pointers to components of state */
  dft_hypers *hyp,		/* Hyperparameters for model */
  int *parents,			/* Parents of nodes for all trees */
  double *divt,			/* Divergence times for all trees */
  double *locations,		/* Node locations for all trees, */
  dft_tree_node *nodes,		/* Internal node representation for all trees */
  int N_cases			/* Number of cases in trees */
)
{
  int dt;

  for (dt = 0; dt<dft->N_trees; dt++)
  { st[dt].d_SD = hyp ? hyp->SD + dft->N_targets*dt : 0;
    st[dt].noise = hyp && model!=0 && model->type=='R' 
                    ? hyp->SD + dft->N_targets*dft->N_trees : 0;
    st[dt].parents = parents ? parents + (2*N_cases)*dt + (N_cases-1) : 0;
    st[dt].divt = divt ? divt + (N_cases-1)*dt - 1 : 0;
    st[dt].locations = locations ? locations+(N_cases-1)*dft->N_targets*dt : 0;
    st[dt].nodes = nodes ? nodes + N_cases*dt : 0;
  }
}


/* SAMPLE HYPERPARAMETER VALUES FROM PRIOR. */

void dft_sample_hyper_prior
( dft_spec *dft,
  model_specification *m,
  dft_hypers *h
)
{
  int dt, j, t;

  j = 0;

  for (dt = 0; dt<dft->N_trees; dt++) 
  { 
    h->d_SD_cm[dt] = prior_pick_sigma (dft->tree[dt].d_SD.width,
                                       dft->tree[dt].d_SD.alpha[0]);
    for (t = 0; t<dft->N_targets; t++)
    { h->SD[j++] = prior_pick_sigma (h->d_SD_cm[dt], 
                                     dft->tree[dt].d_SD.alpha[1]);
    }

    h->c0[dt] = prior_pick_sigma (dft->tree[dt].c0.width,
                                  dft->tree[dt].c0.alpha[0]);
    h->c1[dt] = prior_pick_sigma (dft->tree[dt].c1.width,
                                  dft->tree[dt].c1.alpha[0]);
    h->c2[dt] = prior_pick_sigma (dft->tree[dt].c2.width,
                                  dft->tree[dt].c2.alpha[0]);
  }

  if (m!=0 && m->type=='R')
  {
    h->noise_cm = prior_pick_sigma (m->noise.width, m->noise.alpha[0]);

    for (t = 0; t<dft->N_targets; t++)
    { h->SD[j++] = prior_pick_sigma (h->noise_cm, m->noise.alpha[1]);
    }
  }
}
