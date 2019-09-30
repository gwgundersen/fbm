/* DIST-UTIL.C - Utilities for programs that sample a specified distribution. */

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
#include "mc.h"
#include "formula.h"
#include "data.h"
#include "dist.h"
#include "dist-data.h"


/* COUNT NUMBER OF STATE VARIABLES. */

int dist_count_vars(void)
{
  int i, n;
  char *p;

  n = 0;
  for (p = STATE_VARS; *p; p++)
  { for (i = 0; i<=10; i++)
    { if (formula_var_exists[*p-'a'][i])
      { n += 1;
      }
    }
  }

  return n;
}


/* PACK STATE VARIABLES INTO VECTOR. */

void dist_pack_vars
( double *q
)
{
  char *p;
  int i;

  for (p = STATE_VARS; *p; p++)
  { for (i = 0; i<=10; i++)
    { if (formula_var_exists[*p-'a'][i])
      { *q++ = formula_var[*p-'a'][i]; 
      }
    }
  }
}


/* PACK GRADIENT INTO VECTOR. */

void dist_pack_grad
( double *g
)
{
  char *p;
  int i;

  for (p = STATE_VARS; *p; p++)
  { for (i = 0; i<=10; i++)
    { if (formula_var_exists[*p-'a'][i])
      { *g++ = formula_gradient[*p-'a'][i]; 
      }
    }
  }
}


/* UNPACK STATE VARIABLES FROM VECTOR. */

void dist_unpack_vars
( double *q
)
{
  char *p;
  int i;

  for (p = STATE_VARS; *p; p++)
  { for (i = 0; i<=10; i++)
    { if (formula_var_exists[*p-'a'][i])
      { formula_var[*p-'a'][i] = *q++;
      }
    }
  }
}


/* FIND MINUS LOG OF PRIOR DENSITY FOR BAYESIAN MODEL.  The parameters
   must already have been unpacked into formula_var. */

double dist_prior 
( dist_spec *dst,	/* Specification of the distribution */
  double *grad		/* Place to store gradient, null if not wanted */
)
{ 
  double e;

  e = formula (dst->energy + strlen(dst->energy) + 1,  
               0, 1, grad ? STATE_VARS : 0);

  if (grad) dist_pack_grad(grad);

  return e;
}


/* FIND MINUS THE LOG LIKELIHOOD FOR A TRAINING CASE.  The parameters
   must already have been unpacked into formula_var, and the data must
   already have been read.  The inv_temp parameter is used to multiply
   the likelihood and its gradient; it sholud be set to 1.0 to get the
   usual values. */

double dist_likelihood
( dist_spec *dst,	/* Specification of the distribution */
  int c,		/* Number of training case (from 0) */
  double inv_temp,	/* Factor to multiply log likelihood by */
  double *grad		/* Place to add gradient to, null if not wanted */
)
{ 
  double e;
  int i, t;
  char *p;

  if (c<0 || c>=N_train) abort();

  for (i = 0; i<data_spec->N_inputs && i<10; i++)
  { formula_var['i'-'a'][i] = train_inputs [data_spec->N_inputs*c + i];
  }
  formula_var['i'-'a'][10] = formula_var['i'-'a'][0];

  for (t = 0; t<data_spec->N_targets && t<10; t++)
  { formula_var['t'-'a'][t] = train_targets [data_spec->N_targets*c + t];
  }
  formula_var['t'-'a'][10] = formula_var['t'-'a'][0];

  e = inv_temp * formula (dst->energy, 0, 1, grad ? STATE_VARS : 0);

  if (grad)
  {
    for (p = STATE_VARS; *p; p++)
    { for (i = 0; i<=10; i++)
      { if (formula_var_exists[*p-'a'][i])
        { *grad++ += inv_temp * formula_gradient[*p-'a'][i]; 
        }
      }
    }
  }

  return e;
}


/* FIND MINUS THE TOTAL LOG LIKELIHOOD FOR ALL TRAINING CASES.  The 
   parameters must already have been unpacked into formula_var, and the 
   data must already have been read.  The inv_temp parameter is used to 
   multiply the likelihood and its gradient; it sholud be set to 1.0 to 
   get the usual values. */

double dist_total_likelihood
( dist_spec *dst,	/* Specification of the distribution */
  double inv_temp,	/* Factor to multiply log likelihood by */
  double *grad		/* Place to add gradient to, null if not wanted */
)
{ 
  double e;
  int c;

  e = 0;

  for (c = 0; c<N_train; c++)
  { e += dist_likelihood (dst, c, inv_temp, grad);
  }

  return e;
}


/* SAMPLE FROM THE PRIOR DISTRIBUTION.  Exits with an error if this is not
   possible.  The point drawn is stored in the array passed, and possible
   also in the formula data structure. */

void dist_sample_prior
( dist_spec *dst,	/* Specification of the distribution */
  double *q		/* Place to store values */
)
{ 
  int i, n;

  n = dist_count_vars();

  /* Points from the prior come from standard input if -read-prior set.*/

  if (dst->read_prior)
  {
    for (i = 0; i<n; i++)
    { if (scanf("%lf",&q[i])!=1)
      { if (ferror(stdin))
        { fprintf(stderr,"Error trying to read samples from prior\n");
        }
        else if (feof(stdin))
        { fprintf(stderr,"End of file trying to read samples from prior\n");
        }
        else 
        { fprintf(stderr,"Nonnumeric data trying to read samples from prior\n");
        }
        exit(1);
      }
    }
  }

  /* Otherwise, we try to do it ourselves. */

  else
  { 
    formula_sample (dst->energy+strlen(dst->energy)+1, STATE_VARS);
    dist_pack_vars (q);
  }
}


/* PRINT VALUES ASSOCIATED WITH STATE VARIABLES.  Both the variables and
   the values are found in the formula data structure. */

void dist_print_vars (void)
{ 
  int imax, icnt;
  int lused[26];
  int i, k;
  char *p;

  imax = 0;
  for (p = STATE_VARS; *p; p++)
  { icnt = 0;
    for (i = 0; i<=10; i++)
    { if (formula_var_exists[*p-'a'][i]) 
      { icnt += 1;
      }
    }
    lused[*p-'a'] = icnt>0;
    if (icnt>imax) imax = icnt;
  }

  for (k = 1; k<=imax; k++)
  { for (p = STATE_VARS; *p; p++)
    { if (lused[*p-'a'])
      { icnt = 0;
        for (i = 0; i<=10; i++)
        { if (formula_var_exists[*p-'a'][i])
          { icnt += 1;
            if (icnt==k) break;
          }
        }
        if (icnt==k)
        { if (i==10) printf("    %c = ",*p);
          else       printf("   %c%d = ",*p,i);
          printf("%-10lg",formula_var[*p-'a'][i]);
        }
        else
        { printf("%18s","");
        } 
      }
    }
    printf("\n");
  }
}
