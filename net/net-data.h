/* NET-DATA.H - Interface to module for reading data for networks. */

/* Copyright (c) 1995, 1996 by Radford M. Neal 
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


/* VARIABLES HOLDING TRAINING AND/OR TEST DATA.  When the values or targets
   aren't known, the pointers are null. */

extern data_specifications *data_spec; /* Specifications of data sets */

extern int N_train;		/* Number of training cases */

extern net_values *train_values;/* Values associated with training cases */
extern double *train_targets;	/* True targets for training cases */

extern int N_test;		/* Number of test cases */

extern net_values *test_values;	/* Values associated with test cases */
extern double *test_targets;	/* True targets for test cases */


/* PROCEDURES. */

void net_data_read (int, int, net_arch *, 
                    model_specification *, model_survival *);

void net_data_free (void);
