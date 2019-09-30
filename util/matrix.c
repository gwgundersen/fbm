/* MATRIX.C - Routines for doing matrix operations. */

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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"


/* SET A SQUARE MATRIX TO THE IDENTITY. */

void identity_matrix
( double *m,		/* The matrix to set to the identity */
  int n			/* Dimension of the matrix */
)
{
  int k;

  for (k = n*n-1; k>=0; k--) m[k] = 0;

  for (k = n*n-1; k>=0; k -= n+1) m[k] = 1;
}


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


/* FIND THE TRACE OF THE PRODUCT OF TWO SYMMETRIC MATRICES.  The matrix 
   elements are stored in contiguous locations in row-major order.  Space is 
   present for the full matrices, but only the lower triangles are looked at. */

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
  { s += 2 * inner_product(m1,1,m2,1,i);
    s += m1[i] * m2[i];
    m1 += n;
    m2 += n;
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


/* SOLVE TRIANGULAR SYSTEM USING FORWARD SUBSTITUTION.  Solves Lx=b where
   L is lower-triangular using forward substitution.  The vectors x and b
   may be stored with offsets other than one, as specified. */

void forward_solve
( double *m,		/* The matrix L, only the lower-triangle is looked at */
  double *x,		/* Place to store solution */
  int xo,		/* Offset from one element to the next in x */
  double *b,		/* The vector b */
  int bo,		/* Offset from one element to the next in b */
  int n			/* Dimension of matrix and vectors */
)
{ 
  double *xp;
  int i;

  xp = x;

  for (i = 0; i<n; i++)
  { *xp = (*b - inner_product (m, 1, x, xo, i)) / m[i];
    xp += xo;
    b += bo;
    m += n;
  }
}  


/* SOLVE TRIANGULAR SYSTEM USING BACKWARD SUBSTITUTION.  Solves L'x=b where
   L is lower-triangular using backward substitution.  The vectors x and b
   may be stored with offsets other than one, as specified. */

void backward_solve
( double *m,		/* The matrix L, only the lower-triangle is looked at */
  double *x,		/* Place to store solution */
  int xo,		/* Offset from one element to the next in x */
  double *b,		/* The vector b */
  int bo,		/* Offset from one element to the next in b */
  int n			/* Dimension of matrix and vectors */
)
{ 
  int i;

  x += (n-1)*xo;
  b += (n-1)*bo;
  m += n*n-1;

  for (i = n-1; i>=0; i--)
  { *x = (*b - inner_product (m+n, n, x+xo, xo, (n-1)-i)) / *m;
    x -= xo;
    b -= bo;
    m -= n+1;
  }
}  


/* FIND EIGENVALUES AND EIGENVECTORS BY JACOBI ITERATION.  Finds the 
   eigenvalues and eigenvectors of a symmetric matrix using Jacobi iteration.
   The matrix is passed as the first argument, with only the upper triangle
   being used.  The eigenvalues will appear (in arbitrary order) on the
   diagonal of this matrix, while the off-diagonal elements will be set to
   near zero.  The eigenvectors will be stored as the rows of the second
   matrix passed, unless a null pointer is passed.

   The iteration will proceed until the magnitude of the largest off-diagonal
   element is no greater than the tolerance value passed.  The number of 
   2D rotations done to achieve this is returned as the value of this function.
  
   This procedure tries to keep track of the largest off-diagonal element
   in each row/column.  The largest element according to this record is 
   zeroed in the next rotation.  The record is not entirely accurate, but
   the largest off-diagonal element is guaranteed to be no larger than the
   largest of these records (so the record is safe to use for the termination
   criterion).
*/

int jacobi
( double *m,		/* The matrix, in upper triangle */
  double *e,		/* Place to store eigenvectors, or null */
  double tol,		/* Tolerance ratio */
  int n			/* Size of the matrix */
)
{ 
  double c, s, u, tmp;
  double largest_off;
  double lodi, lodj;
  double mi, mj;
  double ei, ej;
  double x;

  double *lod;
  int *lk;

  int L, i, j, k, l, t;
  int ii, jj, ij;
  int lki, lkj;

  /* Allocate space for recording the largest off-diagonal elements. */

  lod = calloc (n, sizeof *lod);
  lk  = calloc (n, sizeof *lk);

  if (lod==0 || lk==0)
  { fprintf(stderr,"jacobi: Can't allocate space\n");
    exit(1);
  }
 
  /* Initialize the matrix of eigenvectors to the identity. */

  if (e!=0)   
  { identity_matrix (e, n);
  }

  /* Initialize the records of largest off-diagonal elements. */

  for (l = 0; l<n; l++) 
  { lod[l] = 0;
    k = l;
    for (i = 0; i<l; i++)
    { x = m[k];
      if (x>lod[l]) 
      { lod[l] = x;
        lk[l] = k;
      } 
      else if (-x>lod[l]) 
      { lod[l] = -x;
        lk[l] = k;
      } 
      k += n;
    } 
    for (j = l+1; j<n; j++)
    { k += 1;
      x = m[k];
      if (x>lod[l]) 
      { lod[l] = x;
        lk[l] = k;
      } 
      else if (-x>lod[l]) 
      { lod[l] = -x;
        lk[l] = k;
      } 
    } 
  }

  /* Do the iterations, one 2D rotation each time around the loop. */

  for (t = 0; ; t++)
  {
    /* Find the largest off-diagonal element; set i and j to it's location. */

    largest_off = 0;
    L = 0;

    for (l = 1; l<n; l++)
    { if (lod[l]>largest_off)
      { largest_off = lod[l];
        L = l;
      }
    }

    i = lk[L]/n;
    j = lk[L]%n;  
 
    /* Check whether we are done, according to the tolerance criterion. */

    if (largest_off<=tol)
    { free(lod);
      free(lk);
      return t;
    }

    /* Find the locations of the elements undergoing rotation. */
  
    ii = n*i + i;
    jj = n*j + j;
  
    ij = ii + j-i;
  
    mi = m[ii];
    mj = m[jj];
  
    /* Find the sine and cosine of the rotation angle. */
  
    if (m[ij]<1e-100 && m[ij]>-1e-100) /* Set specially to avoid overflow */
    { c = 1;
      s = 0;
    }
    else
    { u = (mj-mi) / (2*m[ij]);
      tmp = u>=0 ? 1/(u+sqrt(1+u*u)) : -1/(-u+sqrt(1+u*u));
      c = 1/sqrt(1+tmp*tmp);
      s = tmp*c;
    }
  
    /* Adjust the matrix being diagonalized by this rotation.  Also keep
       track of the new largest off-diagonal elements for these rows/columns. */
  
    lodi = 0; 
    lodj = 0;
  
    m[ii] = c*c*mi - 2*c*s*m[ij] + s*s*mj;
    m[jj] = s*s*mi + 2*s*c*m[ij] + c*c*mj;
  
    m[ij] = c*s*(mi-mj) + (c*c-s*s)*m[ij];
  
    ii = i;
    jj = j;
  
    for (k = i; k>0; k--)
    { mi = m[ii];
      mj = m[jj];
      m[ii] = c*mi - s*mj;
      if (m[ii]>lodi)       { lodi = m[ii];  lki = ii; } 
      else if (-m[ii]>lodi) { lodi = -m[ii]; lki = ii; }
      m[jj] = c*mj + s*mi;
      if (m[jj]>lodj)       { lodj = m[jj];  lkj = jj; }
      else if (-m[jj]>lodj) { lodj = -m[jj]; lkj = jj; }
      ii += n;
      jj += n;
  
    }
  
    ii += 1;
    jj += n;
  
    for (k = j-i-1; k>0; k--)
    { mi = m[ii];
      mj = m[jj];
      m[ii] = c*mi - s*mj;
      if (m[ii]>lodi)       { lodi = m[ii];  lki = ii; } 
      else if (-m[ii]>lodi) { lodi = -m[ii]; lki = ii; }
      m[jj] = c*mj + s*mi;
      if (m[jj]>lodj)       { lodj = m[jj];  lkj = jj; }
      else if (-m[jj]>lodj) { lodj = -m[jj]; lkj = jj; }
      ii += 1;
      jj += n;
    }
  
    ii += 1;
    jj += 1;
  
    for (k = n-j-1; k>0; k--)
    { mi = m[ii];
      mj = m[jj];
      m[ii] = c*mi - s*mj;
      if (m[ii]>lodi)       { lodi = m[ii];  lki = ii; } 
      else if (-m[ii]>lodi) { lodi = -m[ii]; lki = ii; }
      m[jj] = c*mj + s*mi;
      if (m[jj]>lodj)       { lodj = m[jj];  lkj = jj; }
      else if (-m[jj]>lodj) { lodj = -m[jj]; lkj = jj; }
      ii += 1;
      jj += 1;
    }
  
    lod[i] = lodi;
    lod[j] = lodj;
    lk[i] = lki;
    lk[j] = lkj;
  
    /* Adjust the future eigenvectors by this rotation. */
  
    if (e!=0)
    { 
      ii = n*i;
      jj = n*j;
  
      for (k = n; k>0; k--)
      { ei = e[ii];
        ej = e[jj];
        e[ii] = c*ei - s*ej;
        e[jj] = c*ej + s*ei;
        ii += 1;
        jj += 1;
      }  
    }
  }
}
