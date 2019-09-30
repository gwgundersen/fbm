/* GP-COVM.C - Procedures to compute covariances and their derivatives. */

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

#include "misc.h"
#include "log.h"
#include "data.h"
#include "prior.h"
#include "model.h"
#include "gp.h"
#include "gp-data.h"


/* COMPUTE COVARIANCES BETWEEN SETS OF INPUT POINTS.  Finds the covariances 
   for function values between all the input points in one set (with n1
   points), and all the inputs points in another set (with n2 points).  The
   sets of input points are given by an array in which the inputs for the
   different points follow one another (ie, it's an array in row major order
   in which the rows are points, and the columns are the inputs for those
   points).

   The 'keep' argument should be zero to get the usual covariances.  It may
   be instead set to a pointer to an array of flags indicating which non-linear
   terms of the covariance should be kept; the constant and linear parts are
   never kept if 'keep' is not zero.  This is useful for separating the 
   various components of an additive model.

   The covariances are stored in 'cov', as an array with n1 rows and n2
   columns.  The terms in the covariances for non-linear components may also 
   be stored, in arrays pointed to from exp_cov.  If exp_cov is zero, or if 
   one of the corresponding pointer in exp_covs is zero, the non-linear term 
   will not be stored.  Storing these terms can speed up the computation of the 
   derivatives on the covariance by diff_cov, as well as possibly being useful 
   in themselves.

   Note that the jitter part of the covariance and the noise for a regression 
   model are not included here, since they are not determined by the inputs,
   but rather by the identity of the cases. */

void gp_cov
( gp_spec *gp,		/* Specification for Gaussian process model */
  gp_hypers *h,		/* Values of hyperparameters */
  double *x1,		/* First set of inputs points */
  int n1,		/* Number of points in first set */
  double *x2,		/* Second set of input points */
  int n2,		/* Number of points in second set */
  double *cov,		/* Place to store covariances (n1 by n2 array) */
  double **exp_cov,	/* Places to store exponential terms of covariance;
			   may be zero, or contain zeros, if not desired */
  int *keep		/* Array of flags saying whether to keep non-linear
			   terms in the covariance, or zero to keep them all */
)
{ 
  double r[Max_inputs];
  double d, s, pw, c, v;
  int l, i, j, k1, k2, ni;
  double *p1, *p2, *pc, *ec;
  char *flags;
  int spread;

  ni = gp->N_inputs;

  /* Constant part of covariance. */

  v = gp->has_constant && keep==0 ? exp(2 * *h->constant) : 0;

  pc = cov; 

  for (p1 = x1, k1 = 0; k1<n1; p1 += ni, k1++)
  { for (p2 = x2, k2 = 0; k2<n2; p2 += ni, k2++)
    { *pc++ = v;
    }
  }

  /* Linear part of covariance. */

  if (gp->has_linear && keep==0)
  { 
    flags = gp->linear_flags;
    spread = gp->linear_spread;

    for (i = 0; i<ni; i++)
    { if (!(flags[i]&Flag_omit))
      { r[i] = exp(2 * *h->linear[i]);
        for (j = i-spread; j<=i+spread; j++)
        { if (j!=i && j>=0 && j<ni
           && !(flags[j]&Flag_omit) && (flags[i]&Flag_spread))
          { r[i] += exp(2 * *h->linear[j]);
          }
        }
      }
    }

    pc = cov;

    for (p1 = x1, k1 = 0; k1<n1; p1 += ni, k1++)
    { for (p2 = x2, k2 = 0; k2<n2; p2 += ni, k2++)
      { s = 0;
        for (i = 0; i<ni; i++)
        { if (!(flags[i]&Flag_omit))
          { s += r[i] * p1[i] * p2[i];
          }
        }
        *pc++ += s;
      }
    }
  }

  /* Non-linear parts of covariance. */

  for (l = 0; l<gp->N_exp_parts; l++)
  { 
    if (keep!=0 && !keep[l]) continue;

    flags = gp->exp[l].flags;
    spread = gp->exp[l].spread;

    for (i = 0; i<ni; i++)
    { if (!(flags[i]&Flag_omit))
      { r[i] = exp(*h->exp[l].rel[i]);
      }
    }

    c = 2 * *h->exp[l].scale;
    pw = gp->exp[l].power;

    pc = cov;
    ec = exp_cov==0 ? 0 : exp_cov[l];

    for (p1 = x1, k1 = 0; k1<n1; p1 += ni, k1++)
    { for (p2 = x2, k2 = 0; k2<n2; p2 += ni, k2++)
      { 
        if (pw==2)
        { s = 0;
          for (i = 0; i<ni; i++)
          { if (!(flags[i]&Flag_omit))
            { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
              s += d * d;
            }
          }
        }
        else if (pw==1)
        { s = 0;
          for (i = 0; i<ni; i++)
          { if (!(flags[i]&Flag_omit))
            { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
              if (d<0) d = -d;
              s += d;
            }
          }
        }
        else if (pw<0)
        { if (pw!=-1) abort();
          s = 0;
          for (i = 0; i<ni; i++)
          { if (!(flags[i]&Flag_omit))
            { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
              s += log (1 + d * d);
            }
          }
        }
        else
        { s = 0;
          for (i = 0; i<ni; i++)
          { if (!(flags[i]&Flag_omit))
            { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
              if (d<0) d = -d;
              s += pow(d,pw);
            }
          }
        }

        v = exp(c-s);

        if (ec) *ec++ = v;
        *pc++ += v;
      }
    }
  }
}


/* COMPUTE DERIVATIVE OF COVARIANCE MATRIX.  Computes the derivatives of the
   covariance matrix for a set of input points, with respect to a single 
   hyperparameter that is pointed to by the wrt argument.  The wrt argument 
   should point to someplace pointed to from h.  Note that when hyperparameters
   are linked together, a single hyperparameter will appear at several places 
   in h.

   Only the upper triangular part of the covariance matrix is computed, though
   the arrays include space for the lower part as well.  This procedure returns
   1 if the wrt hyperparameter is one that affects the covariance directly, and 
   0 if it does not affect the covariance directly (in which case the 
   derivatives are also all set to zero). 

   The exp_cov parameter is an array of pointers to covariance matrices
   derived from each exponential term of the covariance.  These may be
   computed by gp_cov.  Any of the pointers may be zero, in which case 
   the terms are re-computed here.

   Note that the jitter part of the covariance and the noise for a regression 
   model are not included here, since they are not determined by the inputs,
   but rather by the identity of the cases. */

int gp_cov_deriv
( gp_spec *gp,		/* Specification for Gaussian process model */
  gp_hypers *h,		/* Values of hyperparameters */
  double **exp_cov,	/* Exponential parts of covariance (zeros if absent) */
  double *wrt,		/* Pointer to hyperparameter we want derivative w.r.t.*/
  double *x,		/* The set of input points (concatenated as one array)*/
  double *deriv,	/* Place to store derivatives of covariances */
  int n			/* Number of input points */
)
{ 
  double r[Max_inputs], s, d, t, c, pw, v;
  int l, i, j, k, k1, k2, ni;
  double *p1, *p2, *pd, *ec;
  int spread, found;
  char *flags;

  ni = gp->N_inputs;

  /* See if it's the constant part. */

  found = 0;
  v = 0;

  if (gp->has_constant)
  { if (wrt==h->constant)
    { v = 2 * exp(2 * *h->constant);
      found = 1;
    }
  }

  /* Set derivative according to constant part, or to zero. */
  
  pd = deriv;
  
  for (p1 = x, k1 = 0; k1<n; p1 += ni, k1++)
  { pd += k1;
    for (p2 = x+k1*ni, k2 = k1; k2<n; p2 += ni, k2++)
    { *pd++ = v;
    }
  }

  /* See if it's in the linear part. */

  if (gp->has_linear)
  {
    flags = gp->linear_flags;
    spread = gp->linear_spread;

    if (wrt==h->linear_cm && gp->linear.alpha[1]==0)
    { t = 2 * exp(2 * *h->linear_cm);
      for (i = 0; i<ni; i++)
      { if (!(flags[i]&Flag_omit))
        { pd = deriv;
          for (p1 = x, k1 = 0; k1<n; p1 += ni, k1++)
          { pd += k1;
            for (p2 = x+k1*ni, k2 = k1; k2<n; p2 += ni, k2++)
            { if (flags[i]&Flag_spread)
              { for (k = i-spread; k<=i+spread; k++)
                { if (k>=0 && k<ni && !(flags[k]&Flag_omit) 
                                   && (flags[k]&Flag_spread))
                  { *pd += p1[k] * p2[k] * t;
                  }
                }
              }
              else
              { *pd += p1[i] * p2[i] * t;
              }
              pd += 1;
            }
          }
        }
      }
      found = 1;
    }
    else
    { for (i = 0; i<ni; i++)
      { if (!(flags[i]&Flag_omit) && wrt==h->linear[i])
        { t = 2 * exp(2 * *h->linear[i]);
          pd = deriv;
          for (p1 = x, k1 = 0; k1<n; p1 += ni, k1++)
          { pd += k1;
            for (p2 = x+k1*ni, k2 = k1; k2<n; p2 += ni, k2++)
            { if (flags[i]&Flag_spread)
              { for (k = i-spread; k<=i+spread; k++)
                { if (k>=0 && k<ni && !(flags[k]&Flag_omit) 
                                   && (flags[k]&Flag_spread))
                  { *pd += p1[k] * p2[k] * t;
                  }
                }
              }
              else
              { *pd += p1[i] * p2[i] * t;
              }
              pd += 1;
            }
          }
          found = 1;
        }
      }
    }
  }

  /* See if it's in one of the exponential parts. */

  for (l = 0; l<gp->N_exp_parts; l++)
  { 
    /* Quickly see if it's not here. */

    for (j = 0; j<ni && wrt!=h->exp[l].rel[j]; j++) ;

    if (wrt!=h->exp[l].scale && j==ni) continue;

    /* Compute a few things for later use. */
 
    flags = gp->exp[l].flags;

    c = 2 * *h->exp[l].scale;
    pw = gp->exp[l].power;

    if (gp->exp[l].relevance.alpha[1]==0)
    { r[0] = exp(*h->exp[l].rel_cm);
      for (i = 1; i<ni; i++)
      { r[i] = r[0];
      }
    }
    else
    { for (i = 0; i<ni; i++)
      { if (!(flags[i]&Flag_omit))
        { r[i] = exp(*h->exp[l].rel[i]);
        }
      }
    }

    /* See if it's a relevance hyperparameter. */

    if (j<ni)
    { 
      if (exp_cov[l]!=0) /* Previously computed covariances are available */
      {
        for (i = 0; i<ni; i++)
        {
          if (wrt==h->exp[l].rel[i])
          {
            pd = deriv;
            ec = exp_cov[l];

            for (p1 = x, k1 = 0; k1<n; p1 += ni, k1++)
            { 
              pd += k1;
              ec += k1;
      
              for (p2 = x+k1*ni, k2 = k1; k2<n; p2 += ni, k2++)
              { 
                d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);

                if (pw==2)
                { t = 2*d*d;
                }
                else if (pw==1)
                { t = d>0 ? d : -d;
                }
                else if (pw<0)
                { if (pw!=-1) abort();
                  t = 2*d*d / (1+d*d);
                }
                else
                { t = pw * pow (d>0 ? d : -d, pw);
                }

                *pd++ -= t * *ec++;
              } 
            }
          }
        }
      }
      else /* Previously computed covariances are not available */
      {
        pd = deriv;

        for (p1 = x, k1 = 0; k1<n; p1 += ni, k1++)
        { 
          pd += k1;
  
          for (p2 = x+k1*ni, k2 = k1; k2<n; p2 += ni, k2++)
          { 
            if (pw==2)
            { s = 0;
              t = 0;
              for (i = 0; i<ni; i++)
              { if (!(flags[i]&Flag_omit))
                { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
                  d = d*d;
                  s += d;
                  if (wrt==h->exp[l].rel[i])
                  { t -= d;
                  }
                }
              }
              t *= 2;
            }
            else if (pw==1)
            { s = 0;
              t = 0;
              for (i = 0; i<ni; i++)
              { if (!(flags[i]&Flag_omit))
                { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
                  if (d<0) d = -d;
                  s += d;
                  if (wrt==h->exp[l].rel[i])
                  { t -= d;
                  }
                }
              }
            }
            else if (pw<0)
            { if (pw!=-1) abort();
              s = 0;
              t = 0;
              for (i = 0; i<ni; i++)
              { if (!(flags[i]&Flag_omit))
                { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
                  s += log(1+d*d);
                  if (wrt==h->exp[l].rel[i])
                  { t -= 2*d*d / (1+d*d);
                  }
                }
              }
            }
            else
            { s = 0;
              t = 0;
              for (i = 0; i<ni; i++)
              { if (!(flags[i]&Flag_omit))
                { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
                  if (d<0) d = -d;
                  d = pow(d,pw);
                  s += d;
                  if (wrt==h->exp[l].rel[i])
                  { t -= d;
                  }
                }
              }
              t *= pw;
            }
  
            *pd++ += t * exp(c-s);
          }
        }
      }
  
      found = 1;
    }

    /* See if it's the scale hyperparameter. */

    if (wrt==h->exp[l].scale)
    { 
      pd = deriv;
      ec = exp_cov[l];

      for (p1 = x, k1 = 0; k1<n; p1 += ni, k1++)
      { 
        pd += k1;
        if (ec!=0) ec += k1;

        for (p2 = x+k1*ni, k2 = k1; k2<n; p2 += ni, k2++)
        { 
          if (ec!=0) /* Previously computed covariances are available */
          { *pd++ += 2 * *ec++;
          }
          else 
          { if (pw==2)
            { s = 0;
              for (i = 0; i<ni; i++)
              { if (!(flags[i]&Flag_omit))
                { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
                  s += d * d;
                }
              }
            }
            else if (pw==1)
            { s = 0;
              for (i = 0; i<ni; i++)
              { if (!(flags[i]&Flag_omit))
                { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
                  if (d<0) d = -d;
                  s += d;
                }
              }
            }
            else if (pw<0)
            { if (pw!=-1) abort();
              s = 0;
              for (i = 0; i<ni; i++)
              { if (!(flags[i]&Flag_omit))
                { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
                  s += log (1 + d * d);
                }
              }
            }
            else
            { s = 0;
              for (i = 0; i<ni; i++)
              { if (!(flags[i]&Flag_omit))
                { d = r[i] * (flags[i]&Flag_delta ? p1[i]!=p2[i] : p1[i]-p2[i]);
                  if (d<0) d = -d;
                  s += pow(d,pw);
                }
              }
            }
            *pd++ += 2 * exp(c-s);
          }
        }
      }

      found = 1;
    }
  }

  return found;
}


/* FIND COVARIANCE MATRIX FOR TRAINING CASES.  Computes the covariance
   matrices for the latent values associated with training cases and/or
   for the target values of training cases.  The training data taken from 
   the standard place (see gp-data.h).  The exp_cov argument can pick up
   the exponential terms of the covariance for later use, if desired. */

void gp_train_cov
( gp_spec *gp,		/* Specification of Gaussian process model */
  model_specification *m, /* Specification of data model */
  gp_hypers *h,		/* Hyperparameters for Gaussian process model */
  int which_output,	/* Which output we are looking at (from 0) */
  double *noise_var,	/* Case-by-case noise variances, or zero */
  double *latent_cov,	/* Place to store covariance of latent values, or zero*/
  double *target_cov,	/* Place to store covariance of targets, or zero */
  double **exp_cov	/* Places to store exponential terms of covariance;
			   may be zero, or contain zeros, if not desired */
)
{
  int i;

  if (gp==0 || h==0 || latent_cov==0 && target_cov==0) abort();

  /* If covariances of latent values are not desired, store them in target_cov,
     where they will later be converted into target covariances. */

  if (latent_cov==0) latent_cov = target_cov;

  /* Compute covariance matrix for the output without added noise, but
     with jitter, and store in latent_cov. */

  gp_cov (gp, h, train_inputs, N_train, train_inputs, N_train,
          latent_cov, exp_cov, 0);

  if (gp->has_jitter)
  { for (i = 0; i<N_train; i++)
    { latent_cov[i*N_train+i] += exp(2 * *h->jitter);
    }
  }

  /* That's all if we don't want target covariances.  Otherwise, check
     that finding covariance of targets is possible. */

  if (target_cov==0) return;

  if (m!=0 && m->type!='R')
  { fprintf(stderr,
    "Finding covariance of targets is not possible - missing latent values?\n");
    exit(1);
  }

  if (m!=0 && m->type=='R' && m->noise.alpha[2]!=0 && noise_var==0)
  { fprintf(stderr,
      "Need case-by-case noise variances, but they're not available\n");
    exit(1);
  }

  /* Copy latent_cov to target_cov if necessary. */

  if (target_cov!=latent_cov)
  { for (i = N_train*N_train - 1; i>=0; i--) 
    { target_cov[i] = latent_cov[i];
    }
  }

  /* Add noise covariance to target_cov, either from case-by-case noise
     variances, or based on the noise hyperparameter. */

  if (m==0)
  { /* Do nothing - same as zero noise variance. */
  }
  else if (m->noise.alpha[2]!=0)
  { if (m->autocorr) abort();
    for (i = 0; i<N_train; i++) 
    { target_cov[i*N_train+i] += noise_var[gp->N_outputs*i+which_output];
    }
  }
  else
  { double n;
    int k;
    n = exp(2 * *h->noise[which_output]);
    for (i = 0; i<N_train; i++) 
    { target_cov[i*N_train+i] += n;
    }
    if (m->autocorr)
    { for (i = 0; i<N_train; i++)
      { for (k = 1; k<=m->n_autocorr; k++) 
        { if (i-k>=0)      target_cov[i*N_train+(i-k)] += n * m->acf[k-1];
          if (i+k<N_train) target_cov[i*N_train+(i+k)] += n * m->acf[k-1];
        }
      }
    }
  }
}
