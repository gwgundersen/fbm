/* MISC.C - Miscellaneous utility procedures. */

/* Copyright (c) 1995 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 *
 * Some of the "range" facilities are adapted from modifications by
 * Carl Edward Rasmussen, 1995.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "misc.h"


/* ADD NUMBERS REPRESENTED BY THEIR LOGARITHMS.  Computes log(exp(a)+exp(b))
   in such a fashion that it works even when a and b have large magnitude. */

double addlogs
( double a,
  double b
)
{ 
  return a>b ? a + log(1+exp(b-a)) : b + log(1+exp(a-b));
}


/* ALLOCATE SPACE AND CHECK FOR ERROR.  Calls 'calloc' to allocate space,
   and then displays an error message and exits if the space couldn't be
   found. */

void *chk_alloc
( unsigned n,		/* Number of elements */
  unsigned size		/* Size of each element */
)
{ 
  void *p;

  p = calloc(n,size);

  if (p==0)
  { fprintf(stderr,"Ran out of memory (while trying to allocate %d bytes)\n",
      n*size);
    exit(1);
  }

  return p;
}


/* PARSE RANGE SPECIFICATION.  Parses a string of the following form:

       [low][:[high]][%modulus|+number]

   'low' and 'high' must be a non-negative integers.  If present, 'modulus' 
   or 'number' must be a positive integer (in fact, 'number' must be greater
   than one).  If 'low' is omitted, it defaults to one (not zero, even though      zero is a legal value).  If neither 'modulus' nor 'number' are specified, 
   a modulus of one is assumed.  If 'high' is omitted but the colon or a 
   modulus or a number is present, 'high' is set to -1 as a signal to the 
   caller.  If just a single number appears, it is interpreted as a value 
   for 'low', and 'high' is set to -2.  Something must appear (i.e. the empty 
   string is not legal).

   The caller provides pointers to places to store 'low', 'high', and 
   'modulus' or 'number'.  The pointer for modulus/number may be null, in 
   which case a modulus/number specification is illegal.  If a 'number'
   is specified, its negation is stored, to distinguish it from a
   'modulus'.

   Incorrect specifications lead to an error message being displayed, and
   the program being terminated.  The checks include verifying that 
   low<=high when 'high' is specified, but the caller may need to check
   this also after in cases where 'high' is left to default. */

void parse_range
( char *str,		/* String to parse */
  int *low,		/* Place to store 'low' */
  int *high,		/* Place to store 'high', or -1 or -2 if absent */
  int *modulus		/* Place to store 'modulus' or 'number' (negated),
                           or a null pointer if these are illegal */
)
{
  char *s;
  char t;

  s = str;

  if (*s==0) goto error;

  /* Look for value for 'low', or let it default. */

  *low = 1;
  
  if (*s>='0' && *s<='9')
  { 
    *low = 0;

    while (*s>='0' && *s<='9') 
    { *low = 10*(*low) + (*s-'0');
      s += 1;
    }
  }

  /* Look for value for 'high', or signal its absence. */

  *high = *s==0 ? -2 : -1;

  if (*s==':')
  { 
    s += 1;

    if (*s!=0 && *s>='0' && *s<='9')
    {
      *high = 0;

      while (*s>='0' && *s<='9') 
      { *high = 10*(*high) + (*s-'0');
        s += 1;
      }

      if (*high<*low)
      { fprintf(stderr,"High end of range is less than low end: %s\n",str);
        exit(1);
      }
    }
  }

  /* Look for value for 'modulus' or 'number' if such is legal, or let 
     modulus default to one. */

  if (modulus!=0)
  { *modulus = 1;
  }

  t = *s;

  if (*s=='%' || *s=='+')
  { 
    if (modulus==0) goto error;

    s += 1;

    *modulus = 0;

    while (*s>='0' && *s<='9') 
    { *modulus = 10*(*modulus) + (*s-'0');
      s += 1;
    }
    
    if (*modulus==0) goto error;

    if (t=='+') 
    { if (*modulus==1) goto error;
      *modulus = -*modulus;
    }
  }

  /* Check for garbage at end. */

  if (*s!=0) goto error;

  return;

  /* Report error. */

error:

  fprintf(stderr,"Bad range specification: %s\n",str);
  exit(1);
}


/* PARSE RANGE GIVEN IN TERMS OF TIMES.  Like the parse_range procedure
   above, except that 'low' and 'high' are floating-point numbers (typically
   representing quantities of compute time), with the default for 'low'
   being zero. */

void parse_time_range
( char *str,		/* String to parse */
  double *low,		/* Place to store 'low' */
  double *high,		/* Place to store 'high', or -1 or -2 if absent */
  int *modulus		/* Place to store 'modulus', or null if illegal */
)
{
  int point;
  char *s;
  char t;

  s = str;

  if (*s==0) goto error;

  /* Look for value for 'low', or let it default. */

  *low = 0.0; 
  point = 0;
  
  while (*s>='0' && *s<='9' || *s=='.') 
  {
    if (*s=='.') 
    { if (point) goto error; 
      point = 1; 
    }
    else 
    { *low = 10*(*low) + (*s-'0');
      point *= 10;
    }

    s += 1;
  }

  if (point) *low /= point;

  /* Look for value for 'high', or signal its absence. */

  *high = *s==0 ? -2.0 : -1.0;

  if (*s==':')
  { 
    s += 1;

    if (*s!=0 && (*s>='0' && *s<='9' || *s=='.'))
    {
      *high = 0.0;
      point = 0;

      while (*s>='0' && *s<='9' || *s=='.' && !point) 
      {
        if (*s=='.')
        { if (point) goto error; 
          point = 1; 
        }
        else 
        { *high = 10*(*high) + (*s-'0');
          point *= 10;
        }

        s += 1;
      }

      if (point) *high /= point;

      if (*high<*low)
      { fprintf(stderr,"High end of time range is less than low end: %s\n",str);
        exit(1);
      }
    }
  }

  /* Look for value for 'modulus' or 'number' if such is legal, or let 
     modulus default to one. */

  if (modulus!=0)
  { *modulus = 1;
  }

  t = *s;

  if (t=='%' || t=='+')
  { 
    if (modulus==0) goto error;

    s++;
    *modulus = 0;

    while (*s>='0' && *s<='9') 
    { *modulus = 10*(*modulus) + (*s-'0');
      s += 1;
    }
    
    if (*modulus==0) goto error;

    if (t=='+') 
    { if (*modulus==1) goto error;
      *modulus = -*modulus;
    }
  }

  /* Check for garbage at end. */

  if (*s!=0) goto error;

  return;

  /* Report error. */

error:

  fprintf(stderr,"Bad time range specification: %s\n",str);
  exit(1);
}
