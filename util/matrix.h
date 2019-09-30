/* MATRIX.H - Interface to matrix routines. */

/* Copyright (c) 1996 by Radford M. Neal 
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

void identity_matrix (double *, int);

double squared_norm  (double *, int, int);
double inner_product (double *, int, double *, int, int);

void matrix_product (double *, double *, double *, int, int, int);
double trace_of_product (double *, double *, int);

int cholesky (double *, int, double *);
int inverse_from_cholesky (double *, double *, double *, int);

void fill_lower_triangle (double *, int);
void fill_upper_triangle (double *, int);

int jacobi (double *, double *, double, int);
