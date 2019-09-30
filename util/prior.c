/* PRIOR.C - Routines dealing with hierarchical priors specifications. */

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

#include "rand.h"
#include "ars.h"
#include "prior.h"


/* PARSE PRIOR SPECIFICATION.  Returns one if successful, zero if not.  
   Stores the values found in the prior structure supplied. */

int prior_parse
( prior_spec *pr,	/* Place to store prior */
  char *s		/* String to parse */
)
{ 
  int i;

  pr->two_point = *s!=0 && s[strlen(s)-1]=='!';

  pr->scale = 0;
  if (*s=='x')
  { pr->scale = 1;
    s += 1;
  }

  if (s[0]=='+' && (s[1]==0 || s[1]==':'))
  { pr->width = 1e10;
  }
  else
  { pr->width = atof(s);
    if (pr->width<=0) return 0;
  }

  for (i = 0; i<Max_alphas; i++)
  {
    s = strchr(s,':')!=0 ? strchr(s,':')+1 : "";

    if (*s==':' || *s==0)
    { pr->alpha[i] = 0;
    }
    else
    { pr->alpha[i] = atof(s);
      if (pr->alpha[i]<=0) return 0;
    }
  }

  if (strchr(s,':')!=0) return 0;
 
  return 1;
}


/* SHOW PRIOR SPECIFICATION.  Stores a syntactic representation of the
   given prior in the string supplied (which had better be big enough to
   hold it).  Returns a pointer to this string. */

char *prior_show
( char *s,		/* Place to store string */
  prior_spec pr		/* Prior */
)
{
  char s1[100], s2[100];
  int i;

  for (i = Max_alphas-1; i>=0 && pr.alpha[i]==0; i--) ;

  strcpy (s2, pr.two_point ? "!" : "");

  for ( ; i>=0; i--)
  {
    if (pr.alpha[i]==0)
    { sprintf(s1,":%s",s2);
    }
    else
    { sprintf(s1,":%.2f%s",pr.alpha[i],s2);
    }

    strcpy(s2,s1);
  }
  
  sprintf (s, "%s%.3f%s", pr.scale ? "x" : " ", pr.width, s2);

  return s;
}


/* COMPUTE WIDTH OF PRIOR ACCOUNTING FOR AUTOMATIC SCALING.  Computes the 
   effective width part of the prior after scaling in accord with the number
   of source units.  (Scaling occurs only if the 'scale' option is set.) */

double prior_width_scaled
( prior_spec *pr,	/* Prior to find scaling for */
  int n			/* Number of source units */
)
{
  double scale, alpha;

  if (!pr->scale)
  { return pr->width;
  } 
  else
  { 
    alpha = pr->alpha[2]!=0 ? pr->alpha[2] : pr->alpha[1];

    if (alpha==0) /* Infinite */
    { scale = n;
    }
    else if (alpha>2)
    { scale = n * (alpha/(alpha-2));
    }
    else if (alpha==2)
    { scale = n < 3 ? n : n * log((double)n);
    }
    else
    { scale = pow ((double)n, 2/alpha);
    }

    return pr->width / sqrt(scale);
  }
}


/* GENERATE SIGMA VALUE FROM GAMMA DISTRIBUTION.  The value returned is
   the square root of the inverse of a precision value picked from a
   Gamma distribution with specified mean and shape. */

double prior_pick_sigma
( double sigma,		/* Square root of inverse of mean precision */
  double alpha		/* Shape parameter */
)
{
  double p;

  p = alpha==0 ? 1 : rand_gamma(alpha/2) / (alpha/2);
      
  return sigma/sqrt(p);
}


/* ADAPTIVE REJECTION SAMPLING FROM CONDITIONAL DISTRIBUTION FOR A SIGMA VALUE.
   Draws a random value from the conditional distribution for a sigma that is 
   defined by its top-down prior and by the sum of the lower-level precision 
   values that it controls, using the Adaptive Rejection Sampling method. */

typedef struct { double w, a, a0, a1, s; } logp_data;

static double logp (double l, double *d, void *vp)
{ logp_data *p = vp;
  double t = exp(l); 
  double v;
  *d = p->a/2 - t*p->a0/(2*p->w) + p->a1*p->s/(2*t);
  v = l*p->a/2 - t*p->a0/(2*p->w) - p->a1*p->s/(2*t);
  return v;
}

double cond_sigma
( double width,		/* Width parameter for top-level prior */
  double alpha0,	/* Alpha for top-level prior */
  double alpha1,	/* Alpha for lower-level prior */
  double sum,		/* Sum of lower-level precisions */
  int n			/* Number of lower-level precision values */
)
{
  logp_data data;

  data.w  = 1 / (width * width);
  data.a  = alpha0 - n*alpha1;
  data.a0 = alpha0;
  data.a1 = alpha1;
  data.s  = sum;

  return exp (-0.5*ars(log(data.w),log(1+1/sqrt(alpha0)),logp,&data));
}
