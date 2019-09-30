/* PRIOR.C - Routines dealing with hierarchical priors specifications. */

/* Copyright (c) 1995, 1996 by Radford M. Neal 
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

#include "rand.h"
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
