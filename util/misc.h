/* MISC.H - Interface to miscellaneious utility procedures. */

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


void *chk_alloc (unsigned, unsigned);	/* 'calloc' with error exit for null */

double addlogs (double, double);	/* Add numbers represented by logs */

void parse_range (char *, int *, int *, int *);	
					/* Parse low:high%modulus range spec */

void parse_time_range (char *, double *, double *, int *);
					/* Parse range with times rather than
					   iteration numbers */
