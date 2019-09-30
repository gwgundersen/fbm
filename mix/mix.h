/* MIX.H - Interface to mixture modeling modules. */

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


/* SPECIFICATION AND PRIORS FOR MIXTURE MODEL.  Specifies the characteristics
   of the model, and the priors for the hyperparameters used in it.  The
   prior mean of the offsets and the concentration parameter are stored 
   in prior_spec structures even though they are presently not variable,
   since this will facilitate future extensions.

   Stored in log files under type 'P'.  Changes may invalidate old log files. */

typedef struct
{
  int N_inputs;			/* Number of input variables, always 0 for now*/
  int N_targets;		/* Number of target variables */

  int N_components;		/* Number of components in mixture, 0=infinity*/

  prior_spec con_prior;		/* Prior for concentration hyperparameter */
  prior_spec SD_prior;		/* Prior for standard deviations */
  prior_spec mean_prior;	/* Prior for means */

  int reserved[10];		/* Reserved for future use */

} mix_spec;


/* HYPERPARAMETERS FOR MIXTURE MODEL.  This structure stores the values
   of the hyperparameters for a mixture model.  Space is allocated for the 
   maximum number of targets, even though this is a bit inefficient, and
   would cause problems with using the standard Markov chain operations.  
   Space is included for possible "noise" standard deviations for Gaussian 
   distributions as specified using model_spec. 

   The concentration parameter is stored here, even though at present it
   is fixed, since it may become variable in a future extension.

   This record is stored in log files under type 'S'.  Changes may invalidate
   old log files. */

typedef struct
{
  double con;			/* Concentration for Dirichlet distribution */
				/*   (stored here in unscaled form)         */

  double SD_cm;			/* Common standard deviation for all targets */
  double SD[Max_targets];	/* SD of offsets for each target */

  double mean[Max_targets];	/* Mean offset for each target */

  double noise_cm;		/* Common SD for Gaussian distributions */
  double noise[Max_targets];	/* Gaussian SD for each target */
  
  int reserved[10];		/* Reserved for future use */

} mix_hypers;


/* PROCEDURES. */

void mix_record_sizes (log_gobbled *);

void mix_check_specs_present (mix_spec *, int, model_specification *);

void mix_hyper_init (mix_spec *, model_specification *, mix_hypers *);

void mix_print_hypers     (mix_spec *, model_specification *, mix_hypers *);
void mix_print_params     (mix_spec *, model_specification *, 
                           int, int *, double *, double *);
void mix_print_indicators (int, short *);

void mix_freq (short *, int, int *, int);
void mix_sort (short *, int, int *, int, double *, double *, int);
