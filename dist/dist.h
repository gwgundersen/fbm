/* DIST.H - Interface for modules that sample from a specified distributio. */

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


#define STATE_VARS "uvwxyz"	/* Letters for names of state variables */


/* SPECIFICATION OF THE DISTRIBUTION. 

   The distribution is specified by an "energy" function, which contains
   references to state variables, or by functions giving the prior and
   the likelihood.  The specification structure contains these formulas (as
   character strings), along with definitions of constants (not involving state
   variables) that can be used in this definition.
 
   This is written to the log file as a record of type 'd' with index -1. */

#define Spec_size 3000	/* Max chars in energy/prior/likelihood specification */

typedef struct
{ 
  char energy[Spec_size];  /* Formula giving the energy, or formulas for the
                              likelihood and prior, followed by any constant 
                              definitions.  Strings are separated by nulls, 
                              with two consecutive nulls at the end. */

  int Bayesian;		/* Is this a Bayesian posterior distribution? */

  int zero_temper;	/* Should tempering be done with reference to zero
                           energy distribution, rather than with reference to
                           the standard normal? */

  int read_prior;	/* Read points from the prior from standard input? */

  int reserved[8];	/* Reserved for future use */

} dist_spec;


/* PROCEDURES. */

int dist_count_vars(void);
void dist_pack_vars(double *), 
     dist_pack_grad(double *), 
     dist_unpack_vars(double *);

double dist_prior (dist_spec *, double *);

double dist_total_likelihood (dist_spec *, double, double *),
       dist_likelihood (dist_spec *, int, double, double *);

void dist_sample_prior (dist_spec *, double *);

void dist_print_vars (void);
