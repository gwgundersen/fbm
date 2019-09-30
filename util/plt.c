/* PLT.C - Skeleton of program to plot quantities obtained from a log file. */

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
#include "quantities.h"


/* SKELETON OF PROGRAM TO PLOT DATA FROM A LOG FILE.  This module is linked
   with application-specific routines to produce a program for writing
   quantities obtained from a log file to standard output, in a form 
   suitable for plotting.  (This program does not actually invoke a plot 
   program, however.) */

enum { PLT, TBL, HIST } program_type = PLT;


/* DISPLAY USAGE MESSAGE AND EXIT.  The plt_usage procedure must be defined 
   by the application specific module. */

void usage(void)
{
  extern void plt_usage(char *);

  plt_usage("plt x-axis y-axis { log-file [ range ] }");

  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{ 
  char *x_axis, *y_axis;
  char **ap, **apa;

  int lindex, hindex, index_mod;

  log_file logf;
  log_gobbled logg;

  int very_first;

  quantities_described qd;
  quantities_held *qh_first, *qh, **qhp;

  double *p_x, *p_y;

  int N_on_x, N_on_y;
  int V_on_x, V_on_y;

  int a, i, v;

  /* Look at arguments other than log files and ranges. */

  if (argc<4) usage();

  x_axis = argv[1];
  y_axis = argv[2];

  ap = argv+3;

  for (apa = ap; *apa!=0 && strcmp(*apa,"/")!=0; apa++) ;

  if (ap==apa) usage();

  if (*apa!=0) apa += 1;

  for (a = 0; quant_app_arguments[a]; a++)
  { (*quant_app_arguments[a])(&apa);
  }

  if (*apa!=0) usage();

  /* Go through all the log files and ranges specified. */

  very_first = 1;

  while (*ap!=0 && strcmp(*ap,"/")!=0)
  {
    /* Look at arguments giving log file and range. */

    logf.file_name = *ap++;

    if (*ap!=0 && strchr(":%0123456789",**ap)!=0)
    { parse_range(*ap,&lindex,&hindex,&index_mod);
      if (index_mod<0)
      { fprintf(stderr,"Bad range specification: %s\n",*ap);
        exit(1);
      }
      ap += 1;
    }
    else
    { lindex = 1;
      hindex = -1;
      index_mod = 1;
    }

    if (hindex<0) hindex = 1000000000;

    /* Open log file and set up for gobbling. */
  
    log_file_open(&logf,0);

    log_gobble_init(&logg,!very_first);
    
    for (a = 0; quant_app_record_sizes[a]; a++)
    { (*quant_app_record_sizes[a])(&logg);
    }

    /* Gobble up records with negative indexes. */

    while (logf.header.index<0)
    { log_gobble(&logf,&logg);
    }

    /* Skip to start of range, gobble up records at start, and let 
       application specific module have a look. */

    while (!logf.at_end 
       && (logf.header.index<lindex || logf.header.index%index_mod!=0))
    { log_file_forward(&logf);
    }

    if (logf.at_end) continue;

    log_gobble(&logf,&logg);

    for (a = 0; quant_app_initialize[a]; a++)
    { (*quant_app_initialize[a])(&logg);
    }

    /* Look at quantities to plot. */

    N_on_x = quantities_requested(qd,x_axis,1);
    N_on_y = quantities_requested(qd,y_axis,0);

    quantities_available(qd,&logg);

    V_on_x = V_on_y = 0;

    for (i = 0; i<N_on_x; i++)
    { V_on_x += qd[i].high - qd[i].low + 1;
    }
    for ( ; i<N_on_y; i++)
    { V_on_y += qd[i].high - qd[i].low + 1;
    }

    if (V_on_x!=1 && V_on_y!=1 && V_on_x!=V_on_y)
    { fprintf(stderr,
  "Number of values on x-axis (%d) not compatible with number on y-axis (%d)\n",
        V_on_x, V_on_y);
      exit(1);
    }

    qh_first = 0;
    qhp = &qh_first;

    /* Go through all the records in the indicated range. */

    for (;;)
    {
      /* Gather values for this iteration. */

      qh = quantities_storage(qd);
      *qhp = qh;
      qh->next = 0;
      qhp = &qh->next;
 
      quantities_evaluate(qd,qh,&logg);

      /* Skip to next desired index, or to end of range. */
  
      while (!logf.at_end && logf.header.index<=hindex 
               && logf.header.index%index_mod!=0)
      { log_file_forward(&logf);
      }

      if (logf.at_end || logf.header.index>hindex)
      { break;
      }

      /* Gobble up records for next index. */

      log_gobble(&logf,&logg);
    }

    /* Plot data collected. */

    for (v = 0; v<V_on_x || v<V_on_y; v++)
    {
      if (!very_first)
      { printf("\n");
      }

      very_first = 0;

      for (qh = qh_first; qh!=0; qh = qh->next)
      { 
        p_x = qh->value[0]      + (V_on_x==1 ? 0 : v);
        p_y = qh->value[N_on_x] + (V_on_y==1 ? 0 : v);
 
        printf("%20.8e %20.8e\n",*p_x,*p_y);
      }
    }

    fflush(stdout);

    /* Re-initialize for next log file. */

    while (qh_first!=0)
    { qh = qh_first;
      qh_first = qh->next;
      free(qh);
    }

    log_file_close(&logf);

    for (a = 0; quant_app_cleanup[a]; a++)
    { (*quant_app_cleanup[a])();
    }
  }

  exit(0);
}
