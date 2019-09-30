/* MODEL.H - Interface to model specification. */

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


/* MODEL FOR TARGET VARIABLES.  Defines the way in which the distribution
   of the target variables in a case is related to the value of the 
   model function for the case, as computed from the inputs (by a neural
   network, for instance).

   Stored in log files under type 'M'.  Changes may invalidate old log files. */

#define Max_autocorr 35	/* Maximum lag of non-zero autocorrelations */

typedef struct
{ 
  int type;		/* Model used for observed data, zero if no model */
			/* 'B' = binary, 'C' = class, 'R' = real,         */
			/* 'V' = survival data                            */

  prior_spec noise;	/* Prior for noise level in model of real targets */

  int autocorr;		/* 0 = no autocorrelations, 1 = autocorr as specified */
  int n_autocorr;	/* Number of specified autocorrelations */
  float acf[Max_autocorr]; /* Autocorrelations at lags up to n_autocorr */

  int reserved[3];	/* Reserved for future use */

} model_specification;


/* CHARACTERISTICS OF MODEL FOR SURVIVAL DATA.  Used only if the model is
   one for survival data.  Records the form of the hazard function, and for 
   a piece-wise linear hazard function, the times of the vertexes (as specified
   in model-spec).

   Stored in log files under type 'V'.  Changes may invalidate old log files. */

#define Max_time_points 100	/* Maximum time points in piece-wise linear   */
                                /*   approximation to hazard function         */
typedef struct
{
  int hazard_type;		/* Type of hazard function:                   */
	                        /*   'C' = constant, 'P' = piece-wise constant*/

  int reserved[4];		/* Reserved for future use                    */

  int log_time;			/* Present times to network in log form?      */

  double time[Max_time_points+1]; /* List of times used in piece-wise constant*/
                                  /*   representation, terminated by zero     */

} model_survival;


/* PROCEDURES. */

int model_targets (model_specification *, int);
