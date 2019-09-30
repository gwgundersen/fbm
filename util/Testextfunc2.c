/* TESTEXTFUNC2.C - Another program for testing external function interface. */

/* The value of Testextfunc2(a1,a2,...,an) is a1*a2*...*an */

#include <stdio.h>
#include <math.h>

#include "extfunc.h"

ext_header ext_head;
double ext_args[Max_ext_args];

main(void)
{
  double v;
  int i, r;

  for (;;)
  { 
    /* Read the next request for function evaluation.  Exit normally
       on EOF (which indicates the calling program has terminated). */

    r = fread (&ext_head, sizeof ext_head, 1, stdin);

    if (r!=1)
    { exit(ferror(stdin));
    }

    /* See if the request is for something we can't do, or is garbage. */

    if (ext_head.want!=Value 
     || ext_head.n_args<0 
     || ext_head.n_args>Max_ext_args) 
    { abort();
    }

    /* Read the argments of the function. */

    r = fread (ext_args, sizeof (double), ext_head.n_args, stdin);

    if (r!=ext_head.n_args) 
    { exit(1);
    }

    /* Compute the value. */

    v = 1;

    for (i = 0; i<ext_head.n_args; i++)
    { v *= ext_args[i]; 
    }

    /* Write the value. */

    fwrite (&v, sizeof (double), 1, stdout);
    fflush (stdout);
  }
}
