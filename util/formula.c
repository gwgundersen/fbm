/* FORMULA.C - Parse and evaluate a mathematical formula. */

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

#include "formula.h"
#include "extfunc.h"
#include "rand.h"

extern double lgamma(double), digamma(double);


/* CONSTANT PI.  Defined here if not in <math.h>. */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* GLOBAL VARIABLES.  Must match declarations in formula.h */

char formula_var_exists[26][11]; /* Whether variables exist */
double formula_var[26][11];	 /* Values of variables */
double formula_gradient[26][11]; /* Derivatives of expression with 
				    respect to certain variables */


/* TABLE OF EXTERNAL FUNCTIONS. */

static struct		/* Table of external functions */
{ 
  char name[Max_ext_name+1];  /* Name of external function */
  FILE *to, *from;	      /* Pipes to write args to and read value from*/

} ext_funcs[Max_ext_funcs];

static int next_ext_func; /* Next free slot in table */


/* LOCAL VARIABLES. */

static char *form;	/* The formula as a character string */
static int create;	/* Create variables that don't exist? */
static int eval;	/* Compute value of formula? */
static char *gradvars;	/* Variables to compute gradient with respect to, */
 			/*   zero if no gradient is to be computed        */

static int n_gv;	/* Number of variables in gradient vector */
static int c_gv[26*11];	/* Letters for gradient variables */
static int i_gv[26*11];	/* Indexes for gradient variables */

static char *nch;	/* Pointer to next character of formula */

static double expr (double [26][11]),
              term (double [26][11]),
              factor (double [26][11]),
              prefactor (double [26][11]),
              var (double [26][11]),
              number (double [26][11]);

static double gaussian  (int, double, double [26][11]),
              expgamma2 (int, double, double [26][11]),
              expgamma  (int, double, double [26][11]);

static void error(void);
static double ext_func(int,double);


/* MACRO TO MOVE TO NEXT NON-SPACE CHARACTER. */

#define next() do { nch += 1; } while (*nch==' ' || *nch=='\n') 


/* SPLIT DEFINITION INTO VARIABLE AND FORMULA.  Splits a variable definition
   into the formula part (returned), and the two components of the variable
   name.  An error message is printed and the program terminated if the
   variable name given is not valid. The formula is not checked for validity,
   except that the program is aborted if it is not of the form var=formula. */

char *formula_def
( char *def,		/* The definition */
  int *c,		/* Set to initial letter of variable (minus 'a')  */
  int *i		/* Set to digit following letter, 10 if no digit  */
)
{
  char *e;

  if (!strchr(def,'=')) abort();

  *c = def[0] - 'a';
  if (DIGIT(def[1]))
  { *i = def[1] - '0';
    e = def+2;
  }
  else
  { *i = 10;
    e = def+1;
  }
  while (*e==' ' || *e=='\n') e += 1;
  if (*e!='=' || !LOWER(def[0]))
  { *strchr(def,'=') = 0;
    fprintf(stderr,"Invalid variable name: %s\n",def);
    exit(1);
  }

  return e+1;
}


/* PARSE FORMULA, MAYBE CREATE VARIABLES, MAYBE COMPUTE VALUE / GRADIENT.  

   Parses the formula given by the string form0.  If eval0 is set, it 
   returns the value of the formula, otherwise it returns garbage.  If
   gradvars0 is non-zero, it also computes the derivatives of the value with 
   respect to those variables whose initial letter is in gradvars0 and which
   are marked as existing at the beginning.  These derivatives are stored
   in formula_gradient.

   Writes an error message and exits if a syntax error is found.  If create0 
   is zero, it also reports an error if an undefined variable is referenced; 
   if create0 is non-zero, such a variable is instead created (with initial 
   value of 0).  */

double formula
( char *form0,		/* The formula, as a character string */
  int create0,		/* Create variables that don't exist? */
  int eval0,		/* Compute value of formula? */
  char *gradvars0	/* Variables to compute gradient with respect to, */
)			/*   zero if no gradient is to be computed        */
{
  double v;
  char *p;
  int i;

  form = form0;
  create = create0;
  eval = eval0 || gradvars0;
  gradvars = gradvars0;

  if (gradvars)
  { n_gv = 0;
    for (p = gradvars; *p; p++)
    { if (!LOWER(*p)) abort();
      for (i = 0; i<=10; i++)
      { if (formula_var_exists[*p-'a'][i])
        { c_gv[n_gv] = *p-'a';
          i_gv[n_gv] = i;
          n_gv += 1;
        }
      }
    }
  }

  nch = form-1;
  next();

  v = expr(formula_gradient);

  if (*nch!=0) error();

  return v;
}


/* RECURSIVE-DESCENT PARSER/EVALUATOR. */

static double expr
( double gr[26][11]
)
{ 
  double gr2[26][11];
  double v, v2;
  int k;

  if (*nch=='+' || *nch=='-')
  { v = 0;
    if (gradvars)
    { for (k = 0; k<n_gv; k++)
      { gr [c_gv[k]] [i_gv[k]] = 0;
      }
    }
  }
  else
  { v = term(gr);
  }

  while (*nch!=0 && *nch!=')' && *nch!=']' && *nch!='}' && *nch!=',')
  { switch (*nch)
    { case '+':
      { next();
        v2 = term(gr2);
        if (eval) 
        { v += v2;
        }
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] += gr2 [c_gv[k]] [i_gv[k]];
          }
        }
        break;
      }
      case '-':
      { next();
        v2 = term(gr2);
        if (eval) 
        { v -= v2;
        }
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] -= gr2 [c_gv[k]] [i_gv[k]];
          }
        }
        break;
      }
      default: 
      { error();
      }
    }
  }

  return v;
}

static double term
( double gr[26][11]
)
{
  double gr2[26][11];
  double v, v2;
  char op;
  int grnz, k;

  v = factor(gr);

  while (*nch!=0 && *nch!='+' && *nch!='-' 
      && *nch!=')' && *nch!=']' && *nch!='}' && *nch!=',')
  { 
    op = *nch;
    if (op=='*' || op=='/') next();
 
    if (v==0 && gradvars)
    { grnz = 0;
      for (k = 0; k<n_gv && !grnz; k++)
      { if (gr [c_gv[k]] [i_gv[k]] != 0)
        { grnz = 1;
        }
      }
    }

    if (v!=0 || gradvars && grnz) 
    { v2 = factor(gr2);
    }
    else
    { double eval_old = eval;
      eval = 0;
      (void) factor(gr2);
      eval = eval_old;
    }

    if (op=='/')
    { if (gradvars && (v!=0 || grnz))
      { for (k = 0; k<n_gv; k++)
        { gr [c_gv[k]] [i_gv[k]] = gr [c_gv[k]] [i_gv[k]] / v2
                                 - v * gr2 [c_gv[k]] [i_gv[k]] / (v2*v2);
        }
      }
      if (eval && v!=0) 
      { v /= v2;
      }
    }
    else 
    { if (gradvars && (v!=0 || grnz))
      { for (k = 0; k<n_gv; k++)
        { gr [c_gv[k]] [i_gv[k]] = gr [c_gv[k]] [i_gv[k]] * v2
                                 + v * gr2 [c_gv[k]] [i_gv[k]];
        }
      }
      if (eval && v!=0) 
      { v *= v2;
      }
    }
  }

  return v;
}

static double factor
( double gr[26][11]
)
{ 
  double gr2[26][11];
  double v, w;  
  int i, k;

  v = prefactor(gr);

  if (*nch=='^')
  { next();
    w = number(gr2);
    i = floor(w);
    if (i!=w || i<=0) error();
    if (eval)
    { w = v;
      v = 1;
      for (k = 0; k<i-1; k++)
      { v *= w;
      }
      if (gradvars)
      { for (k = 0; k<n_gv; k++)
        { gr [c_gv[k]] [i_gv[k]] *= i * v;
        }
      }
      v *= w;
    }
  }

  return v;
}

static double prefactor
( double gr[26][11]
)
{ 
  double v;
  int k;

  switch (*nch)
  {
    case '(':
    { next();
      v = expr(gr);
      if (*nch!=')') error();
      next();
      break;
    }
    case '[':
    { next();
      v = expr(gr);
      if (*nch!=']') error();
      next();
      break;
    }
    case '{':
    { next();
      v = expr(gr);
      if (*nch!='}') error();
      next();
      break;
    }

    case 'A':
    {
      if (nch[0]=='A' && nch[1]=='b' && nch[2]=='s' && !LOWER(nch[3]))
      { nch += 2;
        next();
        v = prefactor(gr);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { if (v==0) 
            { gr [c_gv[k]] [i_gv[k]] = 0;
            }
            else if (v<0)
            { gr [c_gv[k]] [i_gv[k]] = - gr [c_gv[k]] [i_gv[k]];
            }
          }
        }
        if (eval) v = fabs(v);
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'B': goto unknown_function;

    case 'C':
    {
      if (nch[0]=='C' && nch[1]=='o' && nch[2]=='s' && !LOWER(nch[3]))
      { nch += 2;
        next();
        v = prefactor(gr);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] *= - sin(v);
          }
        }
        if (eval) v = cos(v);
      }

      else if (nch[0]=='C' && nch[1]=='o' && nch[2]=='s' && nch[3]=='h'
        && !LOWER(nch[4]))
      { nch += 3;
        next();
        v = prefactor(gr);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] *= sinh(v);
          }
        }
        if (eval) v = cosh(v);
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'D':
    {
      if (nch[0]=='D' && nch[1]=='e' && nch[2]=='l' && nch[3]=='t' 
            && nch[4]=='a' && !LOWER(nch[5]))
      { nch += 4;
        next();
        v = prefactor(gr);
        if (eval) v = v==0 ? 1 : 0;
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] = 0;
          }
        }
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'E':
    {
      if (nch[0]=='E' && nch[1]=='x' && nch[2]=='p' && nch[3]=='G'
                           && nch[4]=='a' && nch[5]=='m' && nch[6]=='m'
                           && nch[7]=='a' && nch[8]=='2' && !DIGIT(nch[9]))
      { nch += 8;
        next();
        v = expgamma2(0,v,gr);
      }

      else if (nch[0]=='E' && nch[1]=='x' && nch[2]=='p' && nch[3]=='G'
                           && nch[4]=='a' && nch[5]=='m' && nch[6]=='m'
                           && nch[7]=='a' && !LOWER(nch[8]))
      { nch += 7;
        next();
        v = expgamma(0,v,gr);
      }

      else if (nch[0]=='E' && nch[1]=='x' && nch[2]=='p' && !LOWER(nch[3]))
      { nch += 2;
        next();
        v = prefactor(gr);
        if (eval) v = exp(v);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] *= v;
          }
        }
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'F':    
    {
      if (nch[0]=='F' && nch[1]=='r' && nch[2]=='a' && nch[3]=='c'
        && !LOWER(nch[4]))
      { nch += 3;
        next();
        v = prefactor(gr);
        if (eval) v = v>=0 ? fmod(v,1.0) : 1 - fmod(-v,1.0);
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'G': 
    {
      if (nch[0]=='G' && nch[1]=='a' && nch[2]=='u' && nch[3]=='s'
                      && nch[4]=='s' && nch[5]=='i' && nch[6]=='a'
                      && nch[7]=='n' && !LOWER(nch[8]))
      { nch += 7;
        next();
        v = gaussian(0,v,gr);
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'H': goto unknown_function;
    case 'I': goto unknown_function;
    case 'J': goto unknown_function;
    case 'K': goto unknown_function;

    case 'L':
    {
      if (nch[0]=='L' && nch[1]=='G' && nch[2]=='a' && nch[3]=='m' 
       && nch[4]=='m' && nch[5]=='a' && !LOWER(nch[6]))
      { nch += 5;
        next();
        v = prefactor(gr);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] *= digamma(v);
          }
        }
        if (eval) v = lgamma(v);
      }

      else if (nch[0]=='L' && nch[1]=='o' && nch[2]=='g' && nch[3]=='S'
       && nch[4]=='u' && nch[5]=='m' && nch[6]=='E' && nch[7]=='x'
       && nch[8]=='p' && !LOWER(nch[9]))
      { 
        char close;

        nch += 8;
        next();

        if (*nch=='(') close = ')';
        else if (*nch=='[') close = ']';
        else if (*nch=='{') close = '}';
        else error();

        next();

        v = expr(gr);

        while (*nch==',')
        { 
          double gr2[26][11];
          double v2, sum, f, f2;

          next();

          v2 = expr(gr2);
          
          if (eval) sum = v>v2 ? v + log(1+exp(v2-v)) : v2 + log(1+exp(v-v2));

          if (gradvars)
          { f = exp(v-sum);
            f2 = exp(v2-sum);
            for (k = 0; k<n_gv; k++)
            { gr [c_gv[k]] [i_gv[k]] = f * gr [c_gv[k]] [i_gv[k]]
                                     + f2 * gr2 [c_gv[k]] [i_gv[k]];
            }
          }

          v = sum;
        }

        if (*nch!=close) error();
        next();
      }

      else if (nch[0]=='L' && nch[1]=='o' && nch[2]=='g' && !LOWER(nch[3]))
      { nch += 2;
        next();
        v = prefactor(gr);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] /= v;
          }
        }
        if (eval) v = log(v);
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'M': goto unknown_function;

    case 'N':
    {
      if (nch[0]=='N' && nch[1]=='o' && nch[2]=='r' && nch[3]=='m'
                      && nch[4]=='a' && nch[5]=='l' && !LOWER(nch[6]))
      { nch += 5;
        next();
        v = gaussian(0,v,gr);
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'O': goto unknown_function;

    case 'P':
    {
      if (nch[0]=='P' && nch[1]=='i' && (nch[2]==0 || !LOWER(nch[2])))
      { nch += 1;
        next();
        if (eval) v = M_PI;
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] = 0;
          }
        }
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'Q': goto unknown_function;
    case 'R': goto unknown_function;

    case 'S':
    {
       if (nch[0]=='S' && nch[1]=='i' && nch[2]=='g' && nch[3]=='n'
        && !LOWER(nch[4]))
      { nch += 3;
        next();
        v = prefactor(gr);
        if (eval) v = v>0 ? 1 : v<0 ? -1 : 0;
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] = 0;
          }
        }
      }

      else if (nch[0]=='S' && nch[1]=='i' && nch[2]=='n' && !LOWER(nch[3]))
      { nch += 2;
        next();
        v = prefactor(gr);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] *= cos(v);
          }
        }
        if (eval) v = sin(v);
      }

      else if (nch[0]=='S' && nch[1]=='i' && nch[2]=='n' && nch[3]=='h'
        && !LOWER(nch[4]))
      { nch += 3;
        next();
        v = prefactor(gr);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] *= cosh(v);
          }
        }
        if (eval) v = sinh(v);
      }

      else if (nch[0]=='S' && nch[1]=='q' && nch[2]=='r' && nch[3]=='t'
        && !LOWER(nch[4]))
      { nch += 3;
        next();
        v = prefactor(gr);
        if (eval) v = sqrt(v);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] /= 2*v;
          }
        }
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'T':
    {
      if (nch[0]=='T' && nch[1]=='a' && nch[2]=='n' && !LOWER(nch[3]))
      { nch += 2;
        next();
        v = prefactor(gr);
        if (eval) v = tan(v);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] *= (1+v*v);
          }
        }
      }

      else if (nch[0]=='T' && nch[1]=='a' && nch[2]=='n' && nch[3]=='h'
        && !LOWER(nch[4]))
      { nch += 3;
        next();
        v = prefactor(gr);
        if (eval) v = tanh(v);
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] *= (1-v*v);
          }
        }
      }

      else if (nch[0]=='T' && nch[1]=='h' && nch[2]=='e' && nch[3]=='t' 
            && nch[4]=='a' && !LOWER(nch[5]))
      { nch += 4;
        next();
        v = prefactor(gr);
        if (eval) v = v>=0 ? 1 : 0;
        if (gradvars)
        { for (k = 0; k<n_gv; k++)
          { gr [c_gv[k]] [i_gv[k]] = 0;
          }
        }
      }

      else
      { goto unknown_function;
      }

      break;
    }

    case 'U': goto unknown_function;
    case 'V': goto unknown_function;
    case 'W': goto unknown_function;
    case 'X': goto unknown_function;
    case 'Y': goto unknown_function;
    case 'Z': goto unknown_function;
  
    case 'a': case 'b': case 'c': case 'd': case 'e': case 'f': case 'g': 
    case 'h': case 'i': case 'j': case 'k': case 'l': case 'm': case 'n': 
    case 'o': case 'p': case 'q': case 'r': case 's': case 't': case 'u': 
    case 'v': case 'w': case 'x': case 'y': case 'z': 
    { 
      v = var(gr);

      if (*nch!='~') break;

      next();

      if (nch[0]=='N' && nch[1]=='o' && nch[2]=='r' && nch[3]=='m'
                      && nch[4]=='a' && nch[5]=='l' && !LOWER(nch[6]))
      { nch += 5;
        next();
        v = gaussian(1,v,gr);
      }

      else if (nch[0]=='G' && nch[1]=='a' && nch[2]=='u' && nch[3]=='s'
                           && nch[4]=='s' && nch[5]=='i' && nch[6]=='a'
                           && nch[7]=='n' && !LOWER(nch[8]))
      { nch += 7;
        next();
        v = gaussian(1,v,gr);
      }

      else if (nch[0]=='E' && nch[1]=='x' && nch[2]=='p' && nch[3]=='G'
                           && nch[4]=='a' && nch[5]=='m' && nch[6]=='m'
                           && nch[7]=='a' && nch[8]=='2' && !DIGIT(nch[9]))
      { nch += 8;
        next();
        v = expgamma2(1,v,gr);
      }

      else if (nch[0]=='E' && nch[1]=='x' && nch[2]=='p' && nch[3]=='G'
                           && nch[4]=='a' && nch[5]=='m' && nch[6]=='m'
                           && nch[7]=='a' && !LOWER(nch[8]))
      { nch += 7;
        next();
        v = expgamma(1,v,gr);
      }

      else if (UPPER(nch[0]))
      { goto unknown_distribution;
      }

      else
      { error();
      }

      break;
    }

    case '0': case '1': case '2': case '3': case '4': 
    case '5': case '6': case '7': case '8': case '9':
    case '.':
    { v = number(gr);
      break;
    }

    default:
    { error();
    }
  }

  return v;

unknown_function:

  return ext_func(0,v);

unknown_distribution:

  return ext_func(1,v);
}

static double var
( double gr[26][11]
)
{
  int c, i, k;

  c = *nch - 'a';

  if (*(nch+1) && DIGIT(*(nch+1))) 
  { nch += 1;
    i = *nch - '0';
  }
  else
  { i = 10;
  }

  if (*(nch+1) && DIGIT(*(nch+1))) 
  { nch += 1;
    error();
  }

  next();

  if (!formula_var_exists[c][i])
  { if (create)
    { formula_var_exists[c][i] = 1;
      formula_var[c][i] = 0;
    }
    else
    { if (i==10)   
      { fprintf (stderr, "Undefined variable: %c\n", c+'a');
      }
      else 
      { fprintf (stderr, "Undefined variable: %c%d\n", c+'a', i);
      }
      exit(1);
    }
  }

  if (gradvars)
  { for (k = 0; k<n_gv; k++)
    { gr [c_gv[k]] [i_gv[k]] = 0;
    }
    if (strchr(gradvars,c+'a'))
    { gr[c][i] = 1; 
    }
  }

  return formula_var[c][i];
}

static double number
( double gr[26][11]
)
{ 
  char n[101];
  int k, dot;
  double v;
  char *p;

  k = 0;
  dot = 0;
  while (*nch && strchr("0123456789.",*nch))
  { if (k>=100) error();
    n[k++] = *nch;
    if (*nch=='.') 
    { if (dot) error();
      dot = 1;
    }
    nch += 1;
  }
  if (*nch=='E' && *(nch+1) && strchr("0123456789+-",*(nch+1)))
  { if (k>=100) error();
    n[k++] = *nch;
    nch += 1;
    if (*nch=='+' || *nch=='-')
    { if (k>=100) error();
      n[k++] = *nch;
      nch += 1;
    }
    if (!strchr("0123456789",*nch)) error();
    while (*nch && strchr("0123456789",*nch)) 
    { if (k>=100) error();
      n[k++] = *nch;
      nch += 1;
    }
  }
  n[k] = 0;

  v = strtod(n,&p);
  if (*p!=0) error();

  while (*nch==' ' || *nch=='\n') nch += 1;

  if (gradvars)
  { for (k = 0; k<n_gv; k++)
    { gr [c_gv[k]] [i_gv[k]] = 0;
    }
  }

  return v;
}


/* HANDLE DENSITY FUNCTIONS FOR VARIOUS DISTRIBUTIONS.  These functions
   are called after the distribution name has been scanned, with *nch
   being the open bracket around the argument list.  They return the
   value of the log density (if eval is set), and set the gradient (if
   gradvars is set).  The function and the tilde syntax are both supported,
   with the value of the variable being passed in for tilde syntax. */

static double gaussian
( int tilde,		/* Does this call use the tilde syntax? */
  double v,		/* Set to value if tilde syntax used */
  double gr[26][11]	/* Gradient, already set if tilde syntax used */
)
{ 
  double gr2[26][11], gr3[26][11];
  double v2, v3;
  char close;
  int k;

  if (*nch=='(') close = ')';
  else if (*nch=='[') close = ']';
  else if (*nch=='{') close = '}';
  else error();

  next();

  if (!tilde)
  { v = expr(gr);
    if (*nch!=',') error();
    next();
  }

  v2 = expr(gr2);
  if (*nch!=',') error();
  next();

  v3 = expr(gr3);
  if (*nch!=close) error();
  next();

  if (gradvars)
  { for (k = 0; k<n_gv; k++)
    { gr [c_gv[k]] [i_gv[k]] = 
        (v-v2)/v3 * (gr [c_gv[k]] [i_gv[k]] - gr2 [c_gv[k]] [i_gv[k]])
         + gr3 [c_gv[k]] [i_gv[k]] * (1/(2*v3) - (v-v2)*(v-v2)/(2*v3*v3));
    }
  }

  if (eval) v = log(2*M_PI*v3)/2 + (v-v2)*(v-v2)/(2*v3);

  return v;
}

static double expgamma2
( int tilde,		/* Does this call use the tilde syntax? */
  double v,		/* Set to value if tilde syntax used */
  double gr[26][11]	/* Gradient, already set if tilde syntax used */
)
{ 
  double gr2[26][11], gr3[26][11];
  double v2, v3;
  char close;
  int k;

  if (*nch=='(') close = ')';
  else if (*nch=='[') close = ']';
  else if (*nch=='{') close = '}';
  else error();

  next();

  if (!tilde)
  { v = expr(gr);
    if (*nch!=',') error();
    next();
  }

  v2 = expr(gr2);
  if (*nch!=',') error();
  next();

  v3 = expr(gr3);
  if (*nch!=close) error();
  next();

  if (gradvars)
  { for (k = 0; k<n_gv; k++)
    { gr [c_gv[k]] [i_gv[k]] = 
        (exp(v)*(v2/2)/v3 - (v2/2)) * gr [c_gv[k]] [i_gv[k]] + 
        ( gr2 [c_gv[k]] [i_gv[k]] == 0 ? 0 
          : (-0.5*log((v2/2)/v3) - 0.5 + digamma(v2/2)/2 + exp(v)/(2*v3) - v/2) 
              * gr2 [c_gv[k]] [i_gv[k]] ) +
        ((v2/2)/v3 - exp(v)*(v2/2)/(v3*v3)) * gr3  [c_gv[k]] [i_gv[k]];
    }
  }

  if (eval) v = - (v2/2)*log((v2/2)/v3) + lgamma(v2/2) 
                + exp(v)*(v2/2)/v3 - (v2/2)*v;

  return v;
}

static double expgamma
( int tilde,		/* Does this call use the tilde syntax? */
  double v,		/* Set to value if tilde syntax used */
  double gr[26][11]	/* Gradient, already set if tilde syntax used */
)
{ 
  double gr2[26][11], gr3[26][11];
  double v2, v3;
  char close;
  int k;

  if (*nch=='(') close = ')';
  else if (*nch=='[') close = ']';
  else if (*nch=='{') close = '}';
  else error();

  next();

  if (!tilde)
  { v = expr(gr);
    if (*nch!=',') error();
    next();
  }

  v2 = expr(gr2);
  if (*nch!=',') error();
  next();

  v3 = expr(gr3);
  if (*nch!=close) error();
  next();

  if (gradvars)
  { for (k = 0; k<n_gv; k++)
    { gr [c_gv[k]] [i_gv[k]] = 
        (exp(v)*v3 - v2) * gr [c_gv[k]] [i_gv[k]] + 
        ( gr2 [c_gv[k]] [i_gv[k]] == 0 ? 0 
           : (- log(v3) + digamma(v2) - v) * gr2 [c_gv[k]] [i_gv[k]] ) +
        (- v2/v3 + exp(v)) * gr3  [c_gv[k]] [i_gv[k]];
    }
  }

  if (eval) v = - v2*log(v3) + lgamma(v2) + exp(v)*v3 - v2*v;

  return v;
}


/* SAMPLE VARIABLES ACCORDING TO FORMULA.  The formula must have the form

      [ + ] v1 ~ D1(...) { + v2 ~ D2(...) }

   where v1, v2, etc. are variables starting with letters in the state
   variables parameter, and D1, D2, etc. are distributions known to this
   module.  State variables that are to be generated must exist beforehand,
   the must be generated exactly once (ie, appear once in the list above),
   and their values cannot be referred to before they are generated. 

   The formula should have been checked for syntax errors already.  It is
   checked here that it is of the form above. */

void formula_sample
( char *form0,		/* Distribution formula */
  char *statevars	/* List of state variables */
)
{
  char sv_exists[26][11];
  double v2, v3;
  int c, i;
  char *p;

  /* Set up so that parts of the formula can be parsed and evaluated. */

  form = form0;
  create = 0;
  eval = 1;
  gradvars = 0;

  /* Remember which state variables exists, but set them to not existing 
     until they are generated. */
 
  for (p = statevars; *p; p++)
  { c = *p-'a';
    for (i = 0; i<=10; i++)
    { sv_exists[c][i] = formula_var_exists[c][i];
      formula_var_exists[c][i] = 0;
    }
  }

  /* Parse formula here down to looking at the expressions that are the
     distribution parameters (which are looked at by expr). */

  nch = form-1;
  next();

  if (*nch=='+') next();

  while (*nch!=0)
  {
    /* Look at variable name */

    if (!LOWER(*nch)) goto form_error;

    c = *nch - 'a';

    if (*(nch+1) && DIGIT(*(nch+1))) 
    { nch += 1;
      i = *nch - '0';
    }
    else
    { i = 10;
    }

    if (*(nch+1) && DIGIT(*(nch+1))) 
    { nch += 1;
      error();
    }

    next();

    if (*nch!='~') goto form_error;

    next();

    if (!strchr(statevars,c+'a'))
    { fprintf(stderr,"Not a state variable: ");
      goto var_error;
    }

    if (!sv_exists[c][i])
    { fprintf(stderr,"Variable does not exist: ");
      goto var_error;
    }

    if (formula_var_exists[c][i])
    { fprintf(stderr,"Variable is generated more than once: ");
      goto var_error;
    }

    if (nch[0]=='N' && nch[1]=='o' && nch[2]=='r' && nch[3]=='m'
                    && nch[4]=='a' && nch[5]=='l' && !LOWER(nch[6])
     || nch[0]=='G' && nch[1]=='a' && nch[2]=='u' && nch[3]=='s'
                    && nch[4]=='s' && nch[5]=='i' && nch[6]=='a'
                    && nch[7]=='n' && !LOWER(nch[8]))
    { 
      nch += nch[0]=='N' ? 5 : 7;
      next();

      if (!strchr("([{",*nch)) abort();
      next();

      v2 = expr(0);
      if (*nch!=',') abort();
      next();

      v3 = expr(0);
      if (!strchr(")]}",*nch)) abort();
      next();

      if (v3<0)
      { fprintf(stderr,
          "Variance for Normal/Gaussian distribution is negative\n");
        exit(1);
      }

      formula_var[c][i] = v2 + sqrt(v3) * rand_gaussian();
    }

    else if (nch[0]=='E' && nch[1]=='x' && nch[2]=='p' && nch[3]=='G'
                         && nch[4]=='a' && nch[5]=='m' && nch[6]=='m'
                         && nch[7]=='a' && nch[8]=='2' && !DIGIT(nch[9]))
    { nch += 8;
      next();

      if (!strchr("([{",*nch)) abort();
      next();

      v2 = expr(0);
      if (*nch!=',') abort();
      next();

      v3 = expr(0);
      if (!strchr(")]}",*nch)) abort();
      next();

      if (v2<0)
      { fprintf(stderr,"Shape parameter for ExpGamma2 is negative\n");
        exit(1);
      }

      if (v3<0)
      { fprintf(stderr,"Mean parameter for ExpGamma2 is negative\n");
        exit(1);
      }

      formula_var[c][i] = log (rand_gamma(v2/2) * v3 / (v2/2));
    }

    else if (nch[0]=='E' && nch[1]=='x' && nch[2]=='p' && nch[3]=='G'
                         && nch[4]=='a' && nch[5]=='m' && nch[6]=='m'
                         && nch[7]=='a' && !LOWER(nch[8]))
    { nch += 7;
      next();

      if (!strchr("([{",*nch)) abort();
      next();

      v2 = expr(0);
      if (*nch!=',') abort();
      next();

      v3 = expr(0);
      if (!strchr(")]}",*nch)) abort();
      next();

      if (v2<0)
      { fprintf(stderr,"Shape parameter for ExpGamma is negative\n");
        exit(1);
      }

      if (v3<0)
      { fprintf(stderr,"Scale parameter for ExpGamma is negative\n");
        exit(1);
      }

      formula_var[c][i] = log (rand_gamma(v2) / v3);
    }

    else if (UPPER(nch[0]))
    { fprintf(stderr,
   "Generation of values from external distributions is not yet implemented\n");
      exit(1);
    }

    else
    { abort();
    }

    formula_var_exists[c][i] = 1;

    if (*nch!=0)
    { if (*nch!='+') goto form_error;
      next();
    }
  }

  /* Check that all state variables have been generated. */

  for (p = statevars; *p; p++)
  { c = *p-'a';
    for (i = 0; i<=10; i++)
    { if (sv_exists[c][i] && !formula_var_exists[c][i])
      { fprintf(stderr,"State variable not generated: ");
        goto var_error;
      }
    }
  }

  return;

  /* Report errors and exit. */

var_error:

  if (i==10)   
  { fprintf (stderr, "%c\n", c+'a');
  }
  else 
  { fprintf (stderr, "%c%d\n", c+'a', i);
  }
  exit(1);

form_error:

  fprintf(stderr,
 "Distribution specification does not have the form v~D(...) + v~D(...) + ...\n"
    );
  exit(1);
}


/* REPORT SYNTAX ERROR. */

static void error(void)
{
  char *p, *s;

  fprintf(stderr,"Syntax error:\n  ");

  p = s = form;

  for (;;)
  { 
    if ((*p==0 || *p=='\n') && p>=nch && s<=nch)
    { fprintf(stderr,"\n  ");
      while (s<nch) 
      { fprintf(stderr,".");
        s += 1;
      }
      fprintf(stderr,"^");
    }

    if (*p==0) break;

    if (*p=='\n')
    { fprintf(stderr,"\n  ");
      s = p+1;
    }
    else
    { fprintf(stderr,"%c",*p);
    }

    p += 1;
  }

  fprintf(stderr,"\n");

# if 0 /* For debugging */
    abort(); 
# endif

  exit(1);
}


/* HANDLE REFERENCE TO EXTERNAL FUNCTION. */

static double ext_func
( int tilde,		/* Is this a reference via tilde notation? */
  double v		/* Set to value if tilde syntax used */
)
{
  char *type = tilde ? "distribution" : "function";

  ext_header ext_head;  
  double ext_args[Max_ext_args]; 
  char name[Max_ext_name+1];
  int pipe_to[2], pipe_from[2];
  char close_bracket;
  int ei, i;

  if (gradvars)
  { fprintf(stderr,"Gradients of external functions aren't implemented yet\n");
    exit(1);
  }

  /* Extract function/distribution name from expression. */

  i = 0;

  while (LOWER(nch[i]) || UPPER(nch[i]) || DIGIT(nch[i]))
  { if (i==Max_ext_name)
    { fprintf(stderr,"External function name is too long (max %d)\n",
              Max_ext_name);
      exit(1);
    }
    name[i] = nch[i];
    i += 1;
  }

  name[i] = 0;

  nch += i-1;
  next();

  /* Note first argument if tilde syntax used. */

  ext_head.n_args = tilde;
  if (tilde)
  { ext_args[0] = v;
  }

  /* Parse and evaluate arguments of function/distribution. */

  if (*nch=='(') close_bracket = ')';
  else if (*nch=='[') close_bracket = ']';
  else if (*nch=='{') close_bracket = '}';
  else error();

  next();

  while (*nch!=close_bracket)
  { 
    if (ext_head.n_args==Max_ext_args)
    { fprintf(stderr,"Too many arguments to external %s (max %d)",
              type, Max_ext_args);
      exit(1);
    }

    v = expr(0);

    ext_args[ext_head.n_args] = v;
    ext_head.n_args += 1;

    if (*nch==',') 
    { next();
    }
    else if (*nch!=close_bracket)
    { error();
    }
  }

  next();

  /* Look for a function by this name in the table. */

  for (ei = next_ext_func-1; ei>=0; ei--)
  { if (strcmp(name,ext_funcs[ei].name)==0) 
    { break;
    }
  }

  /* If function not found in table, add it, start up a new process for it. */

  if (ei<0)
  { if (next_ext_func==Max_ext_funcs)
    { fprintf(stderr,"Too many external functions in use (max %d)\n",
              Max_ext_funcs);
      exit(1);
    }
    ei = next_ext_func;
    next_ext_func += 1;
    strcpy(ext_funcs[ei].name,name);
    if (pipe(pipe_to)!=0 || pipe(pipe_from)!=0)
    { fprintf(stderr,"Can't create pipes for external %s\n",type);
      exit(1);
    }
    switch (fork())
    { case -1: 
      { fprintf(stderr,"Can't create process for external %s\n",type);
        exit(1);
      }
      case 0:			/* child */
      { close(pipe_to[1]);
        close(pipe_from[0]);
        close(0);
        dup(pipe_to[0]);
        close(1);
        dup(pipe_from[1]);
        execl(name,name,(char*)0);
        fprintf(stderr,"Can't run program for external %s %s\n",type,name);
        exit(1);
      }
      default:			/* parent */
      { close(pipe_to[0]);
        close(pipe_from[1]);
        ext_funcs[ei].to = fdopen(pipe_to[1],"wb");
        ext_funcs[ei].from = fdopen(pipe_from[0],"rb");
        break;
      }
    }
  }

  /* Send request to evaluate this function/distribution. */

  if (eval)
  {
    ext_head.want = gradvars ? Value_and_gradient : Value;
  
    if (fwrite (&ext_head, sizeof ext_head, 1, ext_funcs[ei].to) != 1
     || fwrite (ext_args, sizeof (double), ext_head.n_args, 
                ext_funcs[ei].to) != ext_head.n_args
     || fflush (ext_funcs[ei].to) == EOF)
    { fprintf (stderr, "Error writing arguments for external %s %s\n",
               type, ext_funcs[ei].name);
      exit(1);
    }
  
    if (fread (&v, sizeof (double), 1, ext_funcs[ei].from) != 1)
    { fprintf (stderr, "Error reading value for external %s %s\n",
               type, ext_funcs[ei].name);
      exit(1);
    }
  }
  
  /* Return value. */

  return v;
}
