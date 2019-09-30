/* GP.H - Interface to Gaussian Process modules. */

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


/* SPECIFICATION AND PRIORS FOR GAUSSIAN PROCESS MODEL.  Specifies the form 
   of the covariance function, and the priors for the hyperparameters that
   define it. 

   Stored in log files under type 'P'.  Changes may invalidate old log files. */

#define Max_exp_parts 10	/* Max. number of exp parts in the covariance */

#define Flag_omit 1		/* Flag bit indicating input ignored */
#define Flag_delta 2		/* Flag bit indicating delta distance */
#define Flag_spread 4		/* Flag saying relevance should be spread out */

typedef struct
{
  int N_inputs;			/* Number of input variables */
  int N_outputs;		/* Number of output variables */

  int has_constant;		/* Is constant part present? */
  prior_spec constant;		/* Prior for constant part of covariance */

  int has_linear;		/* Is linear part present? */
  prior_spec linear;		/* Prior for linear part of covariance */

  int has_jitter;		/* Is jitter part present? */
  prior_spec jitter;		/* Prior for jitter part of covariance*/

  int N_exp_parts;		/* Number of exponential terms in covariance */

  struct 			/* Exponential terms in covariance function */
  { 
    float power;		  /* Power to raise distance to */

    prior_spec scale;		  /* Prior for scale of this contribution */
    prior_spec relevance;	  /* Prior for relevancies of inputs */

    char flags[Max_inputs];	  /* Flags saying how inputs are treated */

    int spread;			  /* Amount to spread relevance outward by */

    int resvd[9];		  /* Reserved for future use */

  } exp[Max_exp_parts];

  char linear_flags[Max_inputs];/* Flags for linear part */
  int linear_spread;		/* Amount of relevance spread for linear part */

  int reserved[49];		/* Reserved for future use */

} gp_spec;


/* HYPERPARAMETERS FOR GAUSSIAN PROCESS MODEL.  This structure stores
   pointers to the hyperparameters for a Gaussian process model, plus
   hyperparameters pertaining to the data model (if any).  All
   references are indirect, even when there is only one value, in order
   that the whole collection can be stored in a contiguous block, and
   so that hyperparameters that are constrained to be equal can point
   to the same location.  The hyperparameters are stored in logarithmic 
   form (ie, as the logarithms of the hyperparameters in "scale" terms).  

   The contiguous block of hyperparameters is stored in log files under 
   type 'S'.  Changes that alter the contents of this block may invalidate 
   old log files. */

typedef struct
{
  int total_hypers;		/* Total number of hyperparameter values */
  double *hyper_block;		/* Block of all hyperparameter values */

  double *constant;		/* Hyperparameter for constant part */
  double const_constant;	/* Place to store constant value for this part*/

  double *linear_cm;		/* Common hyperparameter for linear part */
  double *linear[Max_inputs];	/* Linear hyperparameters for each input */
  double const_linear;		/* Place to store constant value for this part*/

  double *jitter;		/* Hyperparameter for jitter part */
  double const_jitter;		/* Place to store constant value for this part*/

  struct 			/* Hyperparameters for exponential parts */
  {
    double *scale;		  /* Scale hyperparameter */
    double const_scale;		  /* Place to store constant value for this */

    double *rel_cm;		  /* Common relevance for exp term */
    double *rel[Max_inputs]; 	  /* Relevance hyperparameter for each input */
    double const_rel;		  /* Place to store constant value for this */

  } exp[Max_exp_parts];

  double *noise_cm;		/* Pointer to common sigma for all outputs */
  double *noise[Max_targets];	/* Pointer to sigmas for each target */
  double const_noise;		/* Place to store constant value for this */
  
} gp_hypers;


/* PROCEDURES. */

void gp_record_sizes (log_gobbled *);
void gp_check_specs_present (gp_spec *, int, model_specification *, 
                             model_survival *);

int  gp_hyper_count    (gp_spec *, model_specification *);
void gp_hyper_pointers (gp_hypers *, gp_spec *, model_specification *);

void gp_prior_generate (gp_hypers *, gp_spec *, model_specification *, int,
                        double, double);
double gp_log_prior    (gp_hypers *, gp_spec *, model_specification *, int);
void gp_prior_grad     (gp_hypers *, gp_spec *, model_specification *, 
                        gp_hypers *);

double gp_gdens (double, double, double, int);
void   gp_gdiff (double, double, double, double*, double*);

void gp_cov (gp_spec *, gp_hypers *, double *, int, double *, int, 
             double *, double **, int *);

int gp_cov_deriv (gp_spec *, gp_hypers *, double **, double *, 
                  double *, double *, int);

void gp_train_cov (gp_spec *, model_specification *, gp_hypers *, int,
                   double *, double *, double *, double **);

double gp_likelihood (gp_hypers *, model_specification *, data_specifications *,
                      double *, double *, double *);
