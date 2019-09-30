/* NET-UTIL.C - Various utility procedures for use in neural network code. */

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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "misc.h"
#include "net.h"
#include "rand.h"


/* GENERATE SIGMA VALUE FROM GAMMA DISTRIBUTION.  The value returned is
   the square root of the inverse of a precision value picked from a
   Gamma distribution with specified mean and shape. */

net_sigma net_pick_sigma
( net_sigma sigma,	/* Square root of inverse of mean precision */
  double alpha		/* Shape parameter */
)
{
  double p;

  p = alpha==0 ? 1 : rand_gamma(alpha/2) / (alpha/2);
      
  return sigma/sqrt(p);
}
