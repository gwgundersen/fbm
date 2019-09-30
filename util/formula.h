/* FORMULA.H - Interface to module for handling mathematical formulas. */

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


/* SETS OF CHARACTERS. */

#define LOWER(c) ((c)>='a' && (c)<='z')	/* Lower case letters (for vars)  */
#define UPPER(c) ((c)>='A' && (c)<='Z')	/* Upper case letters (for funcs) */
#define DIGIT(c) ((c)>='0' && (c)<='9')	/* Digits (for modifiers in vars) */


/* TABLES HOLDING INFORMATION ON VARIABLES.  Tables are indexed by the offset
   of the initial letter, and the offset of the following digit, with an offset
   of 10 for a variable consisting of a letter only. */

extern char formula_var_exists[26][11];	/* Whether variables exist */
extern double formula_var[26][11];	/* Values of variables */

extern double formula_gradient[26][11];	/* Derivatives of expression with 
					   respect to certain variables */

/* PROCEDURES. */

double formula (char *, int, int, char *); /* Parse/evaluate formula/gradient */
char *formula_def (char *, int *, int *);  /* Split variable definition */
void formula_sample (char *, char *);      /* Sample variables according to
					      formual of form v~D(...) + ... */
