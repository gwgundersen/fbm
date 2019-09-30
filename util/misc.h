/* MISC.H - Interface to miscellaneious utility procedures. */

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


void *chk_alloc (unsigned, unsigned);	/* 'calloc' with error exit for null */

double addlogs (double, double);	/* Add numbers represented by logs */
double sublogs (double, double);	/* Subtract numbers repr. by logs */

void parse_range (char *, int *, int *, int *);	
					/* Parse low:high%modulus range spec */

void parse_time_range (char *, double *, double *, int *);
					/* Parse range with times rather than
					   iteration numbers */

void parse_flags (char *, char *, int, int);
int list_flags (char *, int, int, char *);
int not_omitted (char *, int, int);
