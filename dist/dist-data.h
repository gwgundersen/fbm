/* DIST-DATA.H - Interface for reading data for Bayesian model. */

/* Copyright (c) 1998 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */


/* VARIABLES HOLDING TRAINING DATA.  The inputs and targets for training
   and test cases are stored with all the inputs or targets for one case
   together - for example, input i for training case c is stored here in 
   train_inputs[N_inputs*c+i]. */

extern data_specifications *data_spec; /* Specifications of data sets */

extern int N_train;		/* Number of training cases */

extern double *train_inputs;	/* Inputs for training cases */
extern double *train_targets;	/* True targets for training cases */


/* PROCEDURES. */

void dist_data_read (void);
void dist_data_free (void);
