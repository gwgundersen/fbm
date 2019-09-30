/* QUANTITIES.C - Routines supporting module for evaluating quantities in log */

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

#include "log.h"
#include "quantities.h"


/* This module handles the interface to application-specific routines for
   extracting the values of various quantities from log files, and defines
   the generic quantities that are available (unless redefined) in all
   applications.  */


/* PARSE STRING TO FIND WHICH QUANTITIES ARE DESIRED.  Parses the quantity
   descriptions in the string passed, and stores the desired quantities in
   the 'qd' structure.  The number of quantity descriptions present in the
   structure is returned.  If the 'init' argument is one, the structure is
   first initialized; otherwise, the quantity descriptions are added to
   those already present.  Errors in the specifications result in a message
   being displayed and the program being terminated. */

int quantities_requested
( quantities_described qd,	/* Place to store which quantities are wanted */
  char *s,			/* String to parse */
  int init			/* Initialize structure first? */
)
{
  int i, n;

  if (init)
  { for (i = 0; i<Max_quantities; i++)
    { qd[i].letter = 0;
      qd[i].available = 0;
      qd[i].modifier = -1;
      qd[i].low = -1;
      qd[i].high = -1;
    }
  }

  for (n = 0; n<Max_quantities && qd[n].letter!=0; n++) ;

  while (*s)
  { 
    if (n==Max_quantities)
    { fprintf(stderr,"Too many quantities specified (max %d)\n",Max_quantities);
      exit(1);
    }

    qd[n].letter = *s++ & 0177;
   
    if (*s>='0' && *s<='9')
    { qd[n].modifier = 0;
      while (*s>='0' && *s<='9')
      { qd[n].modifier *= 10;
        qd[n].modifier += *s - '0';
        s += 1;
      }
    }

    if (*s=='@')
    { 
      s += 1;
      qd[n].low = 0;

      while (*s>='0' && *s<='9')
      { qd[n].low *= 10;
        qd[n].low += *s - '0';
        s += 1;
      }

      qd[n].high = *(s-1)=='@' ? -1 : qd[n].low;

      if (*s==':')
      {
        s += 1;

        qd[n].high = -1;

        if (*s>='0' && *s<='9')
        { 
          qd[n].high = 0;
          while (*s>='0' && *s<='9')
          { qd[n].high *= 10;
            qd[n].high += *s - '0';
            s += 1;
          }

          if (qd[n].high<qd[n].low)
          { fprintf(stderr,"Range for quantity '%c' is bad (%d,%d)\n",
              qd[n].letter, qd[n].low, qd[n].high);
            exit(1);
          }
        }
      }
    }

    n += 1;
  }

  return n;
}


/* FIND OUT WHICH QUANTITIES ARE AVAILABLE.  Writes a message and exits if
   some aren't.  Knows about 't' and '#' itself; otherwise just passes the 
   buck to the application-specific modules.  This procedure also sets upper 
   bound for arrays if that was left indefinite, either as specified by
   the application, or as equal to the low bound. */

void quantities_available 
( quantities_described qd,	/* Descriptions of quantities wanted */
  log_gobbled *logg		/* Records gobbled from log file */
)
{ 
  int a, i;

  for (a = 0; quant_app_available[a]!=0; a++)
  { (*quant_app_available[a])(qd,logg);
  }

  for (i = 0; i<Max_quantities; i++)
  { 
    if (qd[i].letter)
    { 
      if (qd[i].letter=='t' && qd[i].modifier==-1 && qd[i].available==0)
      { qd[i].available = qd[i].low==-1 ? 1 : -1;
      }

      else if (qd[i].letter=='#' && qd[i].available==0)
      { qd[i].available = 
          qd[i].modifier==-1 && qd[i].low!=-1 && qd[i].high!=-1 ? 1 : -1;
      }

      if (qd[i].available<=0)
      { if (qd[i].available==0)
        { fprintf(stderr,"Quantity unknown:");
        }
        else
        { fprintf(stderr,"Quantity not allowed in context:");
        }
        if (qd[i].modifier<0)
        { fprintf(stderr," %c\n", qd[i].letter);
        }
        else 
        { fprintf(stderr," %c%d\n", qd[i].letter, qd[i].modifier);
        }
        exit(1);
      }

      if (qd[i].low!=-1 && qd[i].high==-1)
      { qd[i].high = qd[i].low;
      }

      if (qd[i].high<qd[i].low)
      { fprintf(stderr,"Range for quantity '%c' is bad (%d,%d)\n",
          qd[i].letter, qd[i].low, qd[i].high);
        exit(1);
      }
    }
  }
}


/* ALLOCATE STORAGE FOR QUANTITIES.  The pointer returned is to a single
   allocated block that contains both the quantities_held structure and
   the arrays of doubles that it refers to, which are in order (this allows
   all values to be easily accessed in order). */

quantities_held *quantities_storage
( quantities_described qd	/* Description of quantities being evaluated */
)
{ 
  quantities_held *qh;
  int s, t, i;
  double *b;
  void *v;

  t = 0;

  for (i = 0; i<Max_quantities; i++)
  { if (qd[i].letter)
    { if (qd[i].high<qd[i].low) abort();
      t += qd[i].high - qd[i].low + 1;
    }
  }

  s = sizeof(quantities_held) + t*sizeof(double);
  v = malloc(s);

  if (v==0)
  { fprintf(stderr,"Not enough memory!\n");
    exit(1);
  }

  qh = (quantities_held*) v;
  b  = (double*) (qh+1);

  for (i = 0; i<Max_quantities; i++)
  { if (qd[i].letter)
    { qh->value[i] = b;
      b += qd[i].high - qd[i].low + 1;
      qh->updated[i] = 0;
    }
  }

  return qh;
}


/* EVALUATE QUANTITIES.  Evaluates the 't' quantity (the current index) 
   and the '#' array (with values equal to indexes) itself; otherwise, 
   just passes the buck to the application-specific modules. */

void quantities_evaluate 
( quantities_described qd,	/* Descriptions of quantities desired */
  quantities_held *qh,		/* Storage for values of quantities */
  log_gobbled *logg		/* Gobbled up records to get values from */
)
{
  int a, i, j;

  for (i = 0; i<Max_quantities; i++)
  { qh->updated[i] = 0;
  }

  for (a = 0; quant_app_evaluate[a]!=0; a++)
  { (*quant_app_evaluate[a])(qd,qh,logg);
  }

  for (i = 0; i<Max_quantities; i++)
  { 
    if (qd[i].letter=='t' && qd[i].modifier==-1 && !qh->updated[i])
    { for (j = qd[i].low; j<=qd[i].high; j++)
      { qh->value [i] [j - qd[i].low] = logg->last_index;
      }
      qh->updated[i] = 1;
    }

    else if (qd[i].letter=='#' && !qh->updated[i])
    { for (j = qd[i].low; j<=qd[i].high; j++)
      { qh->value [i] [j-qd[i].low] = j;
      }
      qh->updated[i] = 1;
    }
  }

  for (i = 0; i<Max_quantities; i++)
  { if (qd[i].letter && !qh->updated[i])
    { fprintf(stderr,"Wasn't able to evaluate '%c' at iteration %d\n",
        qd[i].letter, logg->last_index);
      exit(1);
    }
  }
}
