/* BVG.H - Interface for bivariate Gaussian application. */

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


/* SPECIFICATION OF BIVARIATE GAUSSIAN DISTRIBUTION. 

   This is written to the log file as a record of type 'B' with index -1. */

typedef struct
{ 
  double std1, std2;	/* Standard deviations */
  double corr;		/* Correlation */

  int rep;		/* Number of replications */

} bvg_spec;
