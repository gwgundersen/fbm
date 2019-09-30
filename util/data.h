/* DATA.H - Interface to system for specifying data sets. */

/* Copyright (c) 1995 by Radford M. Neal 
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


#define Max_inputs 200		/* Maximum number of "input" values allowed */
#define Max_targets 20		/* Maximum number of "target" values allowed */


/* DESCRIPTION OF A DATA TRANSFORMATION.  This record describes how to 
   transform the raw data from the file to the desired form. */

typedef struct
{ char take_log;		/* Whether to first take the logarithm */
  char data_shift;		/* Is shift amount derived from training data?*/
  char data_scale;		/* Is scale amount derived from training data?*/
  float shift;			/* Amount to add to value after this */
  float scale;			/* Factor to multiply by after shifting */
} data_transformation;


/* RECORD DESCRIBING TRAINING DATA.  This record is stored in the 
   log file, with type 'D' and index -1.  It gives the number of "input"
   and "target" values, the source of the training data in the form of 
   of 'numin' specifications, the transformation to be performed on the raw 
   inputs and targets, and whether the targets are required to be integers, 
   and if so their upper bound (the lower bound is zero). 

   The 'numin' specification for the training and test data should be used in 
   the context of the defauts mentioned in the documentation for 'data-spec.c'
   If not present, the test set specifications are empty strings. */

typedef struct
{
  int N_inputs;			/* Number of "input" values */
  int N_targets;		/* Number of "target" values */
  int int_target;		/* Zero for real targets, else upper bound */

  data_transformation input_trans[Max_inputs];   /* Trans. applied to inputs */
  data_transformation target_trans[Max_targets]; /* Trans. applied to targets */

  char train_inputs[100];	/* Source of training set input values */
  char train_targets[100];	/* Source of training set target values */

  char test_inputs[100];	/* Source of test set input values */
  char test_targets[100];	/* Source of test set target values */

} data_specifications;


/* PROCEDURES. */

data_transformation data_trans_parse (char *);
char *data_trans_build (data_transformation);

double data_trans (double, data_transformation);
double data_inv_trans (double, data_transformation);
