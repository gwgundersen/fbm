/* MOL.H - Interface for molecular dynamics application. */

/* Copyright (c) 1995-2003 by Radford M. Neal 
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


/* SPECIFICATION OF MOLECULAR DYNAMICS SYSTEM.

   This is written to the log file as a record of type 'M' with index -1. */

typedef struct
{ 
  int D;		/* Dimensionality of space (1, 2, or 3) */
  int N;		/* Number of molecules */
  double len_pres;	/* Length of wrapped-around region in each dimension,
                           if positive; minus the pressure if negative */

  double scale;		/* Scale factor for L-J energy */
  double width;		/* Width parameter for L-J energy */
  double max_pair_energy; /* Maximum energy contribution for one pair */
  double inv_temp;	/* Inverse temperature for "zero temperature" 
			   distribution */

  double reserved[10];	/* Reserved for future use */

} mol_spec;


/* MACRO TO WRAP COORDINATES.  Wraps coordinates to the range [0,1), using
   fmod (which returns negative values for negative arguments). */

#define wrap(c) ((c)>=0 ? fmod((c),1.0) : 1.0 - fmod(-(c),1.0))


/* MACRO TO FIND LENGTH OF A DIMENSION.  For the NVT ensemble, this is just 
   a constant.  For the NPT ensemble, the log of the length is stored after
   all the coordinates, but the length is forced to be in the range from
   exp(-100) to exp(100). */

#define dlen(ms,coords) \
  ((ms)->len_pres>0 ? (ms)->len_pres : \
   (coords)[(ms)->N*(ms)->D]>100 ? exp(100) : \
   (coords)[(ms)->N*(ms)->D]<-100 ? exp(-100) : \
   exp((coords)[(ms)->N*(ms)->D]))


/* PROCEDURES. */

double squared_distance  (int, double, double *, int, int);

double energy_and_gradient (double *, double, mol_spec *, int, int, 
                            double *, double *, double *);
