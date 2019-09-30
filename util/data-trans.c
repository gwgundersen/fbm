/* DATA-TRANS.C - Routines for transforming to and from raw data. */

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
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "data.h"


/* PARSE TRANSFORMATION SPECIFICATION.  See data-spec.doc for the syntax
   of a specification.  If the translation and/or scale amounts are to
   be obtained from the training data, the appropriate flags are set,
   but calculation of the amounts is left for the caller. */

data_transformation data_trans_parse
( char *s0
)
{
  data_transformation t;
  char str[100], *s, *p;

  s = s0; 

  t.take_log = 0;
  t.data_shift = 0;
  t.data_scale = 0;
  t.shift = 0;
  t.scale = 1;

  if (strcmp(s,"I")==0) return t;

  if (*s=='L')
  { s += 1;
    t.take_log = 1;
  }

  if (*s=='+' || *s=='-')
  { p = str;
    *p++ = *s++;
    if (*s=='@')
    { t.data_shift = 1;
      s += 1;
    }
    else
    { while (*s!=0 && strchr("0123456789.",*s)!=0) 
      { *p++ = *s++;
      }
      *p = 0;
      t.shift = atof(str);
    }
  }

  if (*s=='x')
  { s += 1;
    if (*s=='@')
    { t.data_scale = 1;
      s += 1;
    }
    else
    { p = str;
      while (*s!=0 && strchr("+-0123456789.",*s)!=0) 
      { *p++ = *s++;
      }
      *p = 0;
      t.scale = atof(str);
    }
  }

  if (*s!=0)
  { fprintf(stderr,"Bad transformation specification: %s\n",s0);
    exit(1);
  }

  return t;
}


/* BUILD DATA TRANSFORMATION SPECIFICATION.  Builds a string specifying   
   the given data transformation.  For translation and scaling derived
   from training data, the actual values are shown, followed by "@".

   The string returned is destroyed by the next call of this procedure. */

char *data_trans_build
( data_transformation t /* Data transformation */
)
{
  static char s[100];
  char u[50];

  s[0] = 0;

  if (t.take_log)
  { strcat(s,"L");
  }

  if (t.data_shift || t.shift!=0)
  { sprintf(u,"%+f",t.shift);
    strcat(s,u);
    if (t.data_shift) strcat(s,"@");
  }

  if (t.data_scale || t.scale!=1)
  { sprintf(u,"x%f",t.scale);
    strcat(s,u);
    if (t.data_scale) strcat(s,"@");
  }

  if (s[0]==0)
  { strcpy(s,"I");
  }

  return s;
}


/* TRANSFORM RAW DATA.  Transforms data from a file to the desired internal
   form.  Reports an error and exits if a non-positive value is to be 
   logarithmically transformed. */

double data_trans
( double d,		/* Data value */
  data_transformation t	/* Transformation to apply */
)
{
  if (isnan(d)) return d;

  if (t.take_log)
  { 
    if (d<=0) 
    { fprintf(stderr,
        "Data to be logarithmically transformed is not positive\n");
      exit(1);
    }

    d = log(d);
  }
  
  return (d+t.shift) * t.scale;
}


/* INVERSE TRANSFORM TO RAW DATA.  Transforms data from its internal form to
   the form in which it appears in the file. */

double data_inv_trans
( double d,		/* Data value */
  data_transformation t	/* Transformation to apply */
)
{
  if (isnan(d)) return d;

  d = (d/t.scale) - t.shift;

  return t.take_log ? exp(d) : d;
}
