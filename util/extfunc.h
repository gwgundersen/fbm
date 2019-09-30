/* EXTFUNC.H - Interface for external functions used in formals. */

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


#define Max_ext_funcs 20 /* Maximum number of external functions in use */
#define Max_ext_name 100 /* Maximum length of name for external function */
#define Max_ext_args 100 /* Maximum number of arguments for external function */

typedef struct		 /* Header written to pipe for external function */
{ 
  enum { Value, Value_and_gradient, Random_variate} want; /* What is wanted */
  int n_args;               /* Number of arguments */

} ext_header;
