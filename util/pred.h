/* PRED.H - Interface between pred program skeleton and application. */

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


/* NAME OF APPLICATION. */

extern char *pred_app_name;	/* Name of application (eg, "net" or "gp") */


/* SOME CONSTANTS. */

#define Median_sample 101	/* Size of median sample per iteration, should
                                   not be too small, or quantiles won't work
				   with few iterations. */

#define Max_median_points 200	/* Max iterations allowed finding median/quant*/


/* SHARED VARIABLES DESCRIBING DATA.  Set by application in pred_app_init. */

extern data_specifications *data_spec; /* Specifications of data sets */

extern int N_test;		/* Number of test cases */

extern double *test_inputs;	/* Inputs for test cases */
extern double *test_targets;	/* True targets for test cases */


/* SHARED VARIABLES SET BY SKELETON MODULE OF PRED PROGRAM.  Looked at by
   the application in pred_app_use_index. */

extern int op_i, op_t, op_r, 	/* Options specified in the xxx-pred command; */
           op_p, op_m, op_n, 	/*   op_N is the number of digit options      */
           op_d, op_l, op_b,    /*   op_R is the 'r' option for count data    */
           op_a, op_q, op_Q,    /*    (for which op_r has been set to 0)      */
           op_D, op_N, op_R;

extern int keep[10];		/* Which components to keep, when op_N != 0 */

extern log_gobbled logg;	/* Records gobbled from log file(s) */

extern model_specification *m;	/* Specification for model */
extern model_survival *sv;	/* Specificaiton for survival model */

extern data_transformation *tr;	/* Transformations applied to inputs/outputs */


extern int M_targets;		/* Number of outputs of model - same as number
				   of targets except for n-way class models */

extern int ms_count;		/* Number of points used to estimate medians */
 
extern int have_targets;	/* Are targets available for test cases? */

extern int use_inverse;		/* Use inverse for computations in gp-pred? */
extern int alt_mean;		/* Use alternate mean computation in gp-pred? */


/* SHARED VARIABLES SET BY THE PRED_APP_USE_INDEX PROCEDURE. */

extern double *test_targ_pred;	/* Mean predictions for test cases at this 
				   iteration; array of size M_targets*N_test */

extern double *test_targ_med;	/* Median predictions for test cases at this 
				   iteration; array of size M_targets*N_test */

extern double *test_log_prob;	/* Log probabilities of test cases based on this
				   iteration - array with N_test entries */

extern float ***median_sample;	/* Holds random sample used to estimate median:
				   array of N_test pointers; each points to an
				   array of N_targets pointers; these point to
				   arrays of Median_sample*Max_median_points */

/* PROCEDURES SUPPLIED BY THE APPLICATION. */

void pred_app_record_sizes (void);
void pred_app_init (void);
int  pred_app_use_index (void);
void pred_app_finish_file (void);
