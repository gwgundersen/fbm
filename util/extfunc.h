/* EXTFUNC.H - Interface for external functions used in formals. */

/* Copyright (c) 2001 by Radford M. Neal 
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


#define Max_ext_funcs 20 /* Maximum number of external functions in use */
#define Max_ext_name 100 /* Maximum length of name for external function */
#define Max_ext_args 100 /* Maximum number of arguments for external function */

typedef struct		 /* Header written to pipe for external function */
{ 
  enum { Value, Value_and_gradient, Random_variate} want; /* What is wanted */
  int n_args;               /* Number of arguments */

} ext_header;
