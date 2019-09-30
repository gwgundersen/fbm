/* MATRIX.C - Routines for doing matrix operations. */

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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"


/* FIND THE SQUARED NORM OF A VECTOR.  The vector (of length n) is
   stored in memory in successive locations that are at a given distance
   apart.  For an ordinary vector, a distance of 1 is used, but other
   distances can come about from looking at columns of a matrix. */

double squared_norm
( double *v,		/* The vector */
  int d,		/* Distance between elements */
  int n			/* Number of elements */
)
{
  double s;
  int i;

  s = 0;

  for (i = n; i>0; i--)
  { s += *v * *v;
    v += d;
  }

  return s;
}


/* FIND THE INNER PRODUCT OF TWO VECTORS.  Each vector is of length n, and
   is stored in memory in successive locations that are at a given distance
   apart.  For an ordinary vector, a distance of 1 is used, but other
   distances can come about from looking at columns of a matrix. */

double inner_product
( double *v1,		/* First vector */
  int d1,		/* Distance between elements of first vector */
  double *v2,		/* Second vector */
  int d2,		/* Distance between elements of second vector */
  int n			/* Number of elements in the vectors */
)
{
  double s;
  int i;

  s = 0;

  for (i = n; i>0; i--)
  { s += *v1 * *v2;
    v1 += d1;
    v2 += d2;
  }

  return s;
}


/* FIND THE PRODUCT OF TWO MATRICES.  The matrix elements are stored in
   contiguous locations in row-major order. */

void matrix_product
( double *m1,		/* Left operand */
  double *m2,		/* Right operand */
  double *r,		/* Place to store result */
  int n,		/* Number of rows in result (and left operand) */
  int m,		/* Number of columns in result (and right operand) */
  int k			/* Number of columns in left operand & rows in right */
)
{
  int i, j;

  for (i = 0; i<n; i++)
  { for (j = 0; j<m; j++)
    { *r++ = inner_product(m1,1,m2+j,m,k);
    }
    m1 += k;
  }
}


/* FIND THE TRACE OF THE PRODUCT OF TWO SQUARE MATRICES.  The matrix elements 
   are stored in contiguous locations in row-major order. */

double trace_of_product
( double *m1,		/* Left operand */
  double *m2,		/* Right operand */
  int n			/* Number of rows & columns in the matrices */
)
{
  double s;
  int i;

  s = 0;

  for (i = 0; i<n; i++)
  { s += inner_product(m1,1,m2,n,n);
    m1 += n;
    m2 += 1;
  }

  return s;
}


/* COMPUTE CHOLESKI DECOMPOSITION, AND DETERMINANT, OF POS-DEF SYMMETRIC MATRIX.
   The Choleski decomosition will overwrite the lower-triangular part of the
   matrix.  The original matrix is assumed to be symmetric, with only the lower
   triangular part actually being looked at.  

   The log of the determinant is stored in the place pointed to by the 
   last parameter, unless this is zero.  The function value is 0 if then
   matrix is not positive definite (or has an eigenvalue very close to zero); 
   the function value is 1 if things were alright.

   The sums used in finding the Cholesky decomposition are computed in the
   reverse of the obvious order.  Because the later values tend to be smaller 
   than the earlier values, this reduces the effect of round-off error. */

int cholesky
( double *m,		/* The matrix */
  int n,		/* Number of rows and columns of matrix */
  double *ld		/* Place to store log of determinant, or zero */
)
{ 
  double *r, *p;
  double s;
  int i, j;

  if (ld) *ld = 0;

  r = m;

  for (i = 0; i<n; i++)
  {
    s = r[i] - squared_norm(r+i-1,-1,i);

    if (s<1e-100)
    { return 0;
    }
    else
    { if (ld) *ld += log(s);
      s = sqrt(s);
      r[i] = s;
      p = r;
      for (j = i+1; j<n; j++)
      { p += n;
        p[i] = (p[i] - inner_product(r+i-1,-1,p+i-1,-1,i)) / s;
      }
    }

    r += n;
  }

  return 1;
}


/* FIND INVERSE FROM CHOLESKY DECOMPOSITION.  Finds the inverse of a 
   symmetric positive definite matrix given its Cholesky decomposition.  
   The inverse will be symmetric, but is nevertheless redundantly stored
   in both the upper and the lower part of the matrix.  The Cholesky
   decomposition that is the input (replaced by the inverse) is in the
   lower triangular part.

   Returns 0 if the cholesky decomposition has a diagonal element less than 
   or equal to zero (or close enough to cause problems); returns 1 if things
   went alright.

   The caller must pass two scratch vectors of length n for temporary
   storage. */

int inverse_from_cholesky
( double *m,		/* Cholesky decomposition, replaced by inverse */
  double *d,		/* First scratch vector */
  double *x,		/* Second scratch vector */
  int n			/* Number of rows and columns in matrix */
)
{
  double *p, *q;
  double s;
  int i, j;

  /* Save the diagonal elements of the Cholesky decomposition in d.  Also check
     that they aren't negative or nearly zero. */

  p = m;
  for (i = 0; i<n; i++) 
  { if (*p<1e-100) 
    { return 0;
    }
    d[i] = *p;
    p += n+1;
  }

  /* Compute inverse one row/column at a time, storing it only in the 
     upper-triangular part of m to avoid destroying the Cholesky 
     decomposition. */

  for (i = 0; i<n; i++)
  { 
    x[i] = 1.0 / d[i];

    p = m + n*i + i;
    for (j = i+1; j<n; j++)
    { p += n;
      x[j] = - inner_product(p,1,x+i,1,j-i) / d[j];
    }

    p = m + i*n + (n-1);
    q = m + n*n + (n-1);
    for (j = n-1; j>=i; j--)
    { s = x[j] - inner_product(q,n,p+1,1,n-j-1);
      *p = s/d[j];
      q -= n+1;
      p -= 1;
    }
  }

  /* Copy the inverse to the symmetric lower part. */

  fill_lower_triangle (m, n);

  return 1;
}


/* FILL IN LOWER TRIANGLE FROM UPPER TRIANGLE.  Sets the elements of
   a square matrix below the diagonal to be the symmetric elements above
   the diagonal. */

void fill_lower_triangle
( double *m,		/* Matrix to fill lower triangle from upper triangle */
  int n			/* Number of rows and columns in matrix */
)
{
  double *p, *q;
  int i, j;

  for (i = 0; i<n; i++)
  { q = m+i;
    p = m+n*i;
    for (j = 0; j<i; j++)
    { *p = *q;
      q += n; 
      p += 1;
    }
  }
}


/* FILL IN UPPER TRIANGLE FROM LOWER TRIANGLE.  Sets the elements of
   a square matrix above the diagonal to be the symmetric elements below
   the diagonal. */

void fill_upper_triangle
( double *m,		/* Matrix to fill upper triangle from lower triangle */
  int n			/* Number of rows and columns in matrix */
)
{
  double *p, *q;
  int i, j;

  for (i = 0; i<n; i++)
  { q = m+i;
    p = m+n*i;
    for (j = 0; j<i; j++)
    { *q = *p;
      q += n; 
      p += 1;
    }
  }
}
