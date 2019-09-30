/* BVG.H - Interface for bivariate Gaussian application. */

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


/* SPECIFICATION OF BIVARIATE GAUSSIAN DISTRIBUTION. 

   This is written to the log file as a record of type 'B' with index -1. */

typedef struct
{ 
  double std1, std2;	/* Standard deviations */
  double corr;		/* Correlation */

  int rep;		/* Number of replications */

} bvg_spec;
