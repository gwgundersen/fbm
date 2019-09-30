/* PRIOR.H - Interface to routines that manipulate hierarchical prior specs. */

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


/* STRUCTURE DESCRIBING HIERARCHICAL PRIOR.  Describes a hierarchical
   prior for a set of hierarchically grouped parameters.  The 'width' field 
   gives the overall scale of the parameters;  it may be scaled according 
   to some quantity (such as the number of units in the source layer of a 
   network).  The parameters are divided into groups (eg, all network weights
   connecting two layers), which may in turn be divided into sub-groups (eg,
   weights on connections out of one unit), and finally each parameters can
   be considered to be a group of one.  The 'alpha' parameters control the 
   hierarchical structure of the prior for these groups, with the alphas
   being the shape parameters of gamma distributions for precisions associated 
   with groups in the hierarchy.  The mean for the top-level precision is
   1/width^2; the upper-level hyperparameters give the mean for the next lower
   layer, etc., until the lowest-level hyperparameters give the precisions
   for Gaussian distributions for the parameters themselves.  Setting an 
   'alpha' to zero means infinity, and effectively eliminates one stage in 
   the hierarchy. 

   Setting the 'two_point' field causes the final distribution to be 
   concentrated at two points (positive and negative), rather than being 
   from a Gaussian.  This is useful for debugging. */

#define Max_alphas 3	/* Maximum depth for hierarchical priors */

typedef struct
{ int scale;		/* Scale width according to, eg, size of source layer?*/
  int two_point;	/* Use two point distribution rather than Gaussian? */
  double width;		/* Value used to determine top-level mean precision */
  double alpha[Max_alphas]; /* Shape parameters down the hierarchy */
} prior_spec;


/* PROCEDURES. */

int prior_parse (prior_spec *, char *);
char *prior_show (char *, prior_spec);

double prior_width_scaled (prior_spec *, int);
double prior_pick_sigma (double, double);

double cond_sigma (double, double, double, double, int);
