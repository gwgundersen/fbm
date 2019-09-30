/* DIST-SPEC.C - Specify distribution to sample from. */

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
#include "dist.h"


static void usage (void), print_formula (char *);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
) 
{
  static dist_spec dspec, *dst = &dspec; /* Initialize to zeros */

  log_file logf;
  log_gobbled logg;

  int c, i, e, j;
  char *a, *f;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing specification. */

  if (argc==2)
  {
    /* Open log file and gobble up initial records. */
  
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    logg.req_size['d'] = sizeof *dst;

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }

    /* Display specification. */  

    if ((dst = logg.data['d'])==0)
    { fprintf(stderr,"No distribution specification in log file\n");
      exit(1);
    }

    printf("\n");

    a = dst->energy + strlen(dst->energy) + 1;
    if (dst->Bayesian)
    { a += strlen(a) + 1;
    }
    if (*a)
    { printf("Constant definitions:\n\n");
      while (*a)
      { print_formula(a);
        a += strlen(a) + 1;
      }
      printf("\n");
    }
    
    if (dst->Bayesian)
    { printf("Minus log of prior density:\n\n");
      print_formula(dst->energy+strlen(dst->energy)+1);
      printf("\n");
      printf("Minus log likelihood:\n\n");
      print_formula(dst->energy);
      printf("\n");
      if (dst->read_prior)
      { printf(
         "Points needed from the prior will be read from standard input\n\n");
      }
    }
    else
    { printf("Energy function:\n\n");
      print_formula(dst->energy);
      printf("\n");
      if (dst->zero_temper)
      { printf("Tempering is with reference to zero energy\n\n"); 
      }
    }

    log_file_close(&logf);

    exit(0);
  }

  /* Otherwise, look at energy/prior/likelihood, and constant definitions */

  dst->Bayesian = 0;
  dst->zero_temper = 0;
  dst->read_prior = 0;

  for (;;)
  { if (argc>1 && strcmp(argv[argc-1],"-zero-temper")==0)
    { dst->zero_temper = 1;
      argc -= 1;
    }
    else if (argc>1 && strcmp(argv[argc-1],"-read-prior")==0)
    { dst->read_prior = 1;
      argc -= 1;
    }
    else
    { break;
    }
  }

  i = 1;
  for (j = 2; j<argc; j++)
  { i += strlen(argv[j])+1;
  }
  if (i>Spec_size)
  { fprintf(stderr,"Distribution specification is too long\n");
    exit(1);
  }

  strcpy(dst->energy,argv[argc-1]);
  a = dst->energy + strlen(dst->energy) + 1;

  if (argc-2>1 && !strchr(argv[argc-2],'='))
  { dst->Bayesian = 1;
    strcpy(a,argv[argc-2]);
    a += strlen(a) + 1;
  }

  if (dst->Bayesian && dst->zero_temper)
  { fprintf(stderr,
"The -zero-temper option cannot be used with Bayesian posterior distributions\n"
    );
    exit(1);
  }

  if (!dst->Bayesian && dst->read_prior)
  { fprintf(stderr,
"The -read_prior option can be used only for Bayesian posterior distributions\n"
    );
    exit(1);
  }

  for (j = 2; j<argc-1-dst->Bayesian; j++)
  { 
    strcpy(a,argv[j]);
    if (!strchr(a,'=')) usage();

    f = formula_def(a,&c,&i);

    if (strchr(STATE_VARS,c+'a'))
    { *strchr(a,'=') = 0;
      fprintf(stderr,"State variable can't be used as a constant: %s\n",a);
      exit(1);
    }

    if (dst->Bayesian && strchr("it",c+'a'))
    { *strchr(a,'=') = 0;
      fprintf(stderr,
        "Input/target variable can't be used as a constant: %s\n",a);
      exit(1);
    }
    
    formula_var[c][i] = formula(f,0,1,0);
    formula_var_exists[c][i] = 2;

    a += strlen(a) + 1;
  }

  *a = 0;

  if (dst->Bayesian) 
  { 
    (void) formula (dst->energy+strlen(dst->energy)+1, 1, 0, 0); /* Prior */

    /* Check that prior refers only constants & state variables. */

    e = 0;  
    for (c = 'a'; c<='z'; c++)
    { if (!strchr(STATE_VARS,c))
      { for (i = 0; i<=10; i++)
        { if (formula_var_exists[c-'a'][i]==1)
          { if (!e) fprintf(stderr,"Invalid state variable:");
            if (i==10)
            { fprintf(stderr," %c",c);
            }
            else
            { fprintf(stderr," %c%d",c,i);
            }
            e = 1;
          }
        }
      }
    }
    if (e) 
    { fprintf(stderr," (state variables must begin with one of %s)\n",
                     STATE_VARS);
      exit(1);
    }
  }

  (void) formula (dst->energy, 1, 0, 0);  /* Log density, or likelihood */

  /* Check that only constants, state variables, and inputs/targets are used. */

  e = 0;  
  for (c = 'a'; c<='z'; c++)
  { if (!strchr(STATE_VARS,c) && (!dst->Bayesian || !strchr("it",c)))
    { for (i = 0; i<=10; i++)
      { if (formula_var_exists[c-'a'][i]==1)
        { if (!e) fprintf(stderr,"Invalid state variable:");
          if (i==10)
          { fprintf(stderr," %c",c);
          }
          else
          { fprintf(stderr," %c%d",c,i);
          }
          e = 1;
        }
      }
    }
  }
  if (e) 
  { fprintf(stderr," (state variables must begin with one of %s)\n",STATE_VARS);
    exit(1);
  }

  /* Create log file and write records. */

  log_file_create(&logf);

  logf.header.type = 'd';
  logf.header.index = -1;
  logf.header.size = sizeof *dst;
  log_file_append(&logf,dst);

  log_file_close(&logf);

  exit(0);
}


/* PRINT FORMULA.  The formula is indented two spaces.  A newline is
   output at the end of the formula.  Two or more spaces are treated
   as equivalent to a newline, since the Bourne shell doesn't cause
   a newline to appear at a continuation using backslash. */

static void print_formula
( char *s
)
{ 
  printf("  ");

  while (*s)
  { if (*s=='\n' || *s==' ' && *(s+1)==' ')
    { printf("\n  ");
      while (*(s+1)==' ') s += 1;
    }
    else
    { printf("%c",*s);
    }
    s += 1;
  }

  printf("\n");
}



/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
"Usage: dist-spec log-file { var=formula } energy [ -zero-temper ]\n");
  fprintf(stderr,
"   or: dist-spec log-file { var=formula } prior likelihood [ -read-prior ]\n");
  fprintf(stderr,
"   or: dist-spec log-file (to display stored specifications)\n");

  exit(1);
}
