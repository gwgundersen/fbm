/* DEMO-QUANTITIES.C - Demo of interface for application specific quantities. */

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
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "log.h"
#include "quantities.h"


/* This is a simple demo of the interface to an application specific module
   for evaluating quantities.  It defines the following quantities:

       e   Has the value 1 if the current index is even, and 0 if odd
       o   Has the value 1 if the current index is odd, and 0 if even
       a   Has a value equal to the application-specific argument passed,
           plus its modifier
*/


/* APPLICATION-SPECIFIC ARGUMENT. */

static char *arg;


/* LOOK AT APPLICATION-SPECIFIC ARGUMETNS.  Passed a pointer to a pointer 
   to an 'argv' style argument list, which should be updated to reflect 
   any argument absorbed. */

void demo_arguments
( char ***argvp
)
{
  if (**argvp!=0)
  {
    arg = **argvp;
    *argvp += 1;
  }
}


/* SPECIFY RECORD SIZES.  Provides the application specific module with an 
   opportunity to say what sizes it expects various records to have. */

void demo_record_sizes
( log_gobbled *logg		/* Records gobbled from log file */
)
{
}


/* INITIALIZE AFTER FIRST RECORDS READ.  Provides the application specific
   module with an opportunity to take an initial look at the records read
   from the log file with indexes less than zero or equal to the first
   index within the range. */

void demo_initialize
( log_gobbled *logg		/* Records gobbled from log file */
)
{
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE.  Looks at 
   the list of quantities and says which ones it will handle, ignoring
   those that are already being handled by another module. */

void demo_available
( quantities_described qd,	/* Description of quantities wanted */
  log_gobbled *logg		/* Records gobbled from log file */
)
{ 
  int i;

  for (i = 0; i<Max_quantities; i++)
  {
    if (qd[i].letter && qd[i].available==0)
    {
      if (qd[i].letter=='e')
      { qd[i].available = qd[i].modifier==-1 && qd[i].low==-1 ? 1 : -1;
      }

      if (qd[i].letter=='o')
      { qd[i].available = qd[i].modifier==-1 && qd[i].low==-1 ? 1 : -1;
      }

      if (qd[i].letter=='a')
      { qd[i].available = arg!=0 && qd[i].modifier>=0 && qd[i].low==-1 ? 1 : -1;
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE.  Evaluates quantities that
   it knows how to evaluate, and which haven't already been evaluated
   by another module. */

void demo_evaluate 
( quantities_described qd, 	/* Description of quantities wanted */
  quantities_held *qh,		/* Place to store values of quantities */
  log_gobbled *logg		/* Records gobbled from log file */
)
{ 
  int i;

  for (i = 0; i<Max_quantities; i++)
  {
    if (qd[i].letter && !qh->updated[i])
    {
      if (qd[i].letter=='e')
      { *qh->value[i] = (logg->last_index+1) % 2;
        qh->updated[i] = 1;
      }
    
      if (qd[i].letter=='o')
      { *qh->value[i] = logg->last_index % 2;
        qh->updated[i] = 1;
      }
    
      if (qd[i].letter=='a')
      { *qh->value[i] = atoi(arg) + qd[i].modifier;
        qh->updated[i] = 1;
      }
    }
  }
}
