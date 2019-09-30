/* DATA.H - Interface to system for specifying data sets. */

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


#define Max_inputs  10000	/* Maximum number of "input" values allowed */
#define Max_targets 10000	/* Maximum number of "target" values allowed */


/* DESCRIPTION OF A DATA TRANSFORMATION.  This record describes how to 
   transform the raw data from the file to the desired form. */

typedef struct
{ char take_log;		/* Whether to first take the logarithm */
  char data_shift;		/* Is shift amount derived from training data?*/
  char data_scale;		/* Is scale amount derived from training data?*/
  float shift;			/* Amount to add to value after this */
  float scale;			/* Factor to multiply by after shifting */
} data_transformation;


/* RECORD DESCRIBING TRAINING AND TEST DATA.  This record is stored in the 
   log file, with type 'D' and index -1.  It gives the number of "input"
   and "target" values, the source of the training data in the form of 
   of 'numin' specifications, the transformations to be performed on the raw 
   inputs and targets, and whether the targets are required to be integers, 
   and if so their upper bound (the lower bound is zero). 

   The 'numin' specification for the training and test data should be used in 
   the context of the defauts mentioned in the documentation for 'data-spec.c'
   If not present, the test set specifications are empty strings. 

   This record is of variable length, varying with the number of inputs and
   targets, since the space needed for transformations varies.  The proper
   size if given by the data_spec_size macro. */

typedef struct
{
  int N_inputs;			/* Number of "input" values */
  int N_targets;		/* Number of "target" values */
  int int_target;		/* Zero for real targets, else upper bound */

  char train_inputs[100];	/* Source of training set input values */
  char train_targets[100];	/* Source of training set target values */

  char test_inputs[100];	/* Source of test set input values */
  char test_targets[100];	/* Source of test set target values */

  data_transformation trans[1];	/* Transformations for inputs and targets, with
				   tr. for inputs first, followed by targets */
} data_specifications;

#define data_spec_size(n_inputs,n_targets) \
( sizeof (data_specifications) \
   + sizeof (data_transformation) * ((n_inputs) + (n_targets) - 1) )


/* PROCEDURES. */

data_transformation data_trans_parse (char *);
char *data_trans_build (data_transformation);

double data_trans (double, data_transformation);
double data_inv_trans (double, data_transformation);
