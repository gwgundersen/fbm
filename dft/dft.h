/* DFT.H - Interface to Dirichlet diffusion tree modules. */

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


/* Specifications for diffusion tree models are stored in log files in a
   record of type 'P', in the format of the dft_spec type below.  This
   specification may be supplemented with one written by model-spec.

   The hyperparameters and parameters at each iteraton (whose number doesn't 
   depend on the number of cases) are stored in a record of type 'S', with 
   the dft_hypers type.

   The tree(s) for an iteration are stored in records of several types, not
   all of which need always be present.  All these records are arrays indexed
   by the number of the tree and the index of a case or a node.  The convention
   used is that trees are numbered from 1, cases from 1, and nodes from -1 
   down for internal nodes, and 1 up for nodes associated with training cases 
   (node 0 does not exist).  The record types and contents (arrays of doubles 
   or ints) are as follows (N is the number of training cases):

      'T'   Divergence times (doubles) of internal nodes for the first tree, 
            followed by divergence times for the second tree, etc.  Times for 
            one tree are for node -1, node -2, ..., up to node -(N-1).

      'R'   Indexes of parents (ints from -1 to -(N-1)) of each node of the
            first tree, followed by the same for the second tree, etc.  The 
            order is parent of node -(N-1), parent of node -(N-2), ..., parent
            of node -1, an empty slot, parent of case 1, parent of case 2, ...
            parent of case N.  The empty slot where the parent of nonexistent
            node 0 would go simplifies indexing.  In the programs, the array 
            for one tree is sometimes references via a pointer to this empty
            slot.

      'L'   Latent vectors for each training case.  Each vector consists of
            as many doubles as there are target variables.  May be absent 
            in some situations.

      'N'   Locations of non-terminal nodes in the first tree, followed by the 
            same for the second tree, etc.  Each location is a vector of
            as many doubles as there are target variables.  The first vector
            is for node -(N-1), then for node -(N-2), ..., up to node -1.
            Locations for terminal nodes are not stored, as they can easily
            be integrated over.  May be absent in some situations.
*/            


/* SPECIFICATION AND PRIORS FOR DIFFUSION TREE MODEL.  Specifies the 
   characteristics of the model, and the priors for parameters used in it.  

   Stored in log files under type 'P'.  Changes may invalidate old log files. */

#define Max_trees 9		/* Maximum number of diffusion trees in model */

typedef struct
{
  int N_inputs;			/* Number of input variables, always 0 for now*/
  int N_targets;		/* Number of target variables */
  int N_trees;			/* Number of diffusion trees in model */

  struct			/* Priors for each tree */
  { 
    prior_spec d_SD;		  /* Prior for diffusion standard deviation */
    prior_spec c0, c1, c2;	  /* Priors for divergence function parameters*/
    int reserved[5];		  /* Reserved for future use */

  } tree[Max_trees];

  int reserved[10];		/* Reserved for future use */

} dft_spec;


/* PARAMETERS FOR DIFFUSION TREE MODEL.  This structure stores the values
   of the overall hyperparameters and of the parameters for each tree.  Space 
   is also included for possible "noise" standard deviations for real-valued
   data, as specified using model_spec.  

   The size of this record varies with the number of targets.

   This record is stored in log files under type 'S'.  Changes may invalidate
   old log files. */

typedef struct
{
  double c0[Max_trees]; 	/* Parameters of divergence functions */
  double c1[Max_trees];
  double c2[Max_trees];

  double d_SD_cm[Max_trees];	/* Common diffusion SD for all targets */
  double noise_cm;		/* Common noise SD for all targets */

  int reserved;			/* Reserved for future use */

  double SD[1];			/* Diffusion SD for each target variable, for
				   each tree, followed by noise SD for each 
				   target, if data is real */
} dft_hypers;


/* REPRESENTATION OF TREE.  Trees are represented in log files using
   indexes of parents, but using indexes of children is more convenient for
   manipulation.  This representation takes the form of an array of nodes 
   of the type below (non-terminal nodes only).  Non-terminal nodes are
   numbered -1, -2, etc.; terminal nodes are numbered 1, 2, etc.  Arrays
   of these nodes are indexed by the absolute value of the node number
   (with index 0 unused). 

   Outside the routines for manipulating these structures, the fields should 
   be accessed only via the chld, npts, and totpts macros defined here, which 
   enforce a consistent ordering of children, regardless of the order in the 
   child array. */

typedef struct
{
  int child[2];			/* The two children, -ve for non-terminals,
				   +ve for terminal nodes (ie, data points) */

  int n_points[2];		/* Number of data points ultimately reached
				   via each child. */
} dft_tree_node;

#define chld(n,i) ((n).child[0]<(n).child[1] ? (n).child[i] \
                                             : (n).child[1-i])

#define npts(n,i) ((n).child[0]<(n).child[1] ? (n).n_points[i] \
                                             : (n).n_points[1-i])

#define totpts(n) ((n).n_points[0] + (n).n_points[1])


/* POINTERS TO STATE OF EACH TREE.  This array of structures holds pointers to 
   the various components of the state pertaining to each tree, offset for 
   easy indexing.  It is indexed by tree numbers that start at 0.  For 
   convenience, a pointer to the noise standard deviations is included as
   well, even though these are not specific to any one tree. */

typedef struct
{
  double *d_SD;		/* Standard dev. for each target var., indexed from 0 */
  double *noise;	/* Noise standard deviations for each target variable */
  int *parents;		/* Parent indexes, offset for use with +/- indexes */
  double *divt;		/* Divergence times, offset for indexing from 1 */
  double *locations;	/* Locations of nodes in tree, or null */
  dft_tree_node *nodes;	/* Nodes in tree, indexed from 1 */

} dft_state[Max_trees];


/* STRUCTURE RECORDING LIKELIHOOD FOR A NODE IN A TREE.  Each likelihood
   function is Gausian in form, and is represented by a mean and variance
   defining its shape, and the log of the value at the peak.  The peak 
   value excludes infinite factors resulting from variance; the number of
   such factors is instead counted in infp. */

typedef struct
{
  double mean;		/* Mean of Gaussian shape of likelihood */
  double var;		/* Variance for Gaussian shape of likelihood */
  double peak;		/* Log of likelihood at peak of Gaussian */
  int infp;		/* Number of infinite density points */

} dft_likelihood;


/* FLOATING-POINT & MATHEMATICAL CONSTANTS. */

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

#define log_sqrt_2pi 0.918938533204672741780329736405


/* PROCEDURES. */

void dft_record_sizes (log_gobbled *);
void dft_check_specs_present (dft_spec *, int, model_specification *);

int  dft_hyper_size  (dft_spec *, model_specification *);
int  dft_divt_size   (dft_spec *, int);
int  dft_parents_size(dft_spec *, int);
int  dft_latent_size (dft_spec *, int);
int  dft_locations_size (dft_spec *, int);
void dft_check_sizes (log_gobbled *, dft_spec *, model_specification *, int);

void dft_hyper_init  (dft_spec *, model_specification *, dft_hypers *);
void dft_setup_state (dft_spec *, model_specification *, dft_state, 
                      dft_hypers *, int *, double *, double *, dft_tree_node *, 
                      int);

int dft_conv_tree   (int *, dft_tree_node *, int);
int dft_add_node    (int *, dft_tree_node *, int, int, int);
int dft_remove_node (int *, dft_tree_node *, int, int);
int dft_sibling     (int *, dft_tree_node *, int);
int dft_max_depth   (int *, int);

void dft_print_hypers (dft_spec *, model_specification *, dft_hypers *);
void dft_print_params (dft_spec *, dft_hypers *);
void dft_print_latent (dft_spec *, double *, int, int);
void dft_print_trees  (dft_spec *, dft_state, int, int);
void dft_graph_trees  (dft_spec *, dft_state, int, double, int);
void dft_print_nodes  (dft_spec *, dft_state, int, int, int);

double dft_div   (dft_hypers *, int, double);
double dft_cdiv  (dft_hypers *, int, double);
double dft_icdiv (dft_hypers *, int, double);

void dft_gen_path (dft_hypers *, dft_state, int, int, int *, double *, int *);

double dft_log_prob_div (dft_spec *, dft_hypers *, int, dft_state, int, double);
double dft_log_prob_paths (dft_spec *, dft_hypers *, int, dft_state, int);
double dft_log_prob_path  (dft_spec *, dft_hypers *, int, dft_state, int);
void   dft_log_prob_node  (dft_spec *, int, dft_state, double *, int, 
                           double, double *, int *, dft_likelihood *);

void dft_node_likelihood      (dft_spec *, int, dft_state, double *, int, int, 
                               double, double *, double *);
void dft_terminal_likelihood  (dft_spec *, int, dft_state, double *, int, int, 
                               double, double *, double *);
void dft_tree_likelihood      (dft_spec *, int, dft_state, double *, int, int, 
                               double, double *, double *, double *, int *,
                               dft_likelihood *);
void dft_update_likelihoods   (dft_spec *, int, dft_state, double *, int,
                               dft_likelihood *);
void dft_combine_likelihoods  (int, double, double, double, int,
                               int, double, double, double, int,
                               int, double *, double, 
                               double *, double *, double *, int *);

void dft_sample_hyper_prior (dft_spec *, model_specification *, dft_hypers *);
void dft_sample_prior     (dft_spec *, dft_state, double *);
void dft_sample_latent    (dft_spec *, model_specification *, dft_state, double,
                           double *);
void dft_gibbs_latent     (dft_spec *, model_specification *, dft_state, double,
                           double *);
void dft_gibbs_locations  (dft_spec *, model_specification *, dft_state, double,
                           double *);
void dft_sample_terminals (dft_spec *, dft_state, double *, 
                           double *[Max_trees]);
