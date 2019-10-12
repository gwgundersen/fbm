/* SRC.H - Interface for source location application. */

/* Copyright (c) 2007 by Radford M. Neal 
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


/* SPECIFICATION OF PRIOR FOR NUMBER AND LOCATION OF SOURCES.

   This is written to the log file as a record of type 'S' with index -1. */

typedef struct
{ 
  int lowN, highN;		/* Range for number of sources */
  double lowQ, highQ;		/* Range for source intensities */
  double powQ;			/* Power transformation for intensities */

  double max_start;		/* Maximum start time for sources */
  double max_stop;		/* Maximum stop time for sources, can be 1e30 */
  double max_duration;		/* Maximum duration for sources, can be 1e30 */

  double low[3];		/* Low end of range for coordinates */
  double high[3];		/* High end of range for coodinates */

  double reserved[20];		/* Reserved for future use */

} src_spec;


/* DETECTOR SPECIFICATION.

   This is written to the log file as a record of type 'T' with index -1. */

typedef struct
{ 
  double log_low_width;		/* Low end of prior range for log of width */
  double log_high_width;	/* High end of prior range for log of width */

  double inv_low_df;		/* Reciprocal of low end of range for df */
  double inv_high_df;		/* Reciprocal of high end of range for df */

  double reserved[20];		/* Reserved for future use */

} det_spec;


/* FLOW SPECIFICATION.

   Currently only two closely related types of flow specification are
   implemented, identified by type 't' (for 'test') and 'T' (for 
   'test-start-stop).

   This is written to the log file as a record of type 'F' with index -1. */

#define steady_state "t"	/* Types that are steady-state models */

typedef struct
{ 
  char type;			/* Type of flow specification */
  char res[3];			/* Reserved for future use */

  double lowU, highU;		/* Prior range for wind speed */

  double ay, by;		/* Parameters in sigma_y formula */
  double az, bz;		/* Parameters in sigma_z formula */

  double reserved[20];		/* Reserved for future use */

} flow_spec;


/* MARKOV CHAIN STATE / MODEL PARAMETERS.  Note that changes here may require
   changes to mc_app_set_range in src-mc.c.

   This is written to the log file as a record of type 'q'. */

typedef struct
{
  double N0;			/* Number of sources (real form)       */

  double U;			/* Wind speed (test flow model)        */

  double log_width;		/* Log of width of noise distribution  */
  double inv_df;		/* Reciprocal of df for noise dist.    */

  struct src_desc		/* Parameters for each source:         */
  { 
    double coord[3];		/*   - x, y, z coordinates             */
    double Q;			/*   - intensity raised to given power */
    double start;		/*   - start time                      */
    double stop;		/*   - stop time                       */

  } src[1]; /* Actual number determined at run-time */

} src_params;


/* PROCEDURES */

int  src_params_length (int);
void src_record_sizes (log_gobbled *);
void src_default_parameters (src_spec *, det_spec *, flow_spec *, src_params *);
void src_check_specs_present (src_spec *, det_spec *, flow_spec *);

void src_prior_generate (src_params *, src_spec *, det_spec *, flow_spec *);
void src_prior_one_src (src_params *, src_spec *, det_spec *, flow_spec *, int);

double src_total (src_spec *, flow_spec *, src_params *, double [4]);
double src_test_cstar (flow_spec *, src_params *, double [3], double [3]);
double src_testss_cstar 
  (flow_spec *, src_params *, double [3], double, double, double [4]);
