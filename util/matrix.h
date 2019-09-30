/* MATRIX.H - Interface to matrix routines. */

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

void identity_matrix (double *, int);

double squared_norm  (double *, int, int);
double inner_product (double *, int, double *, int, int);

void matrix_product (double *, double *, double *, int, int, int);
double trace_of_product (double *, double *, int);

int cholesky (double *, int, double *);
int inverse_from_cholesky (double *, double *, double *, int);

void fill_lower_triangle (double *, int);
void fill_upper_triangle (double *, int);

void forward_solve (double *, double *, int, double *, int, int);
void backward_solve (double *, double *, int, double *, int, int);

int jacobi (double *, double *, double, int);
