/* HIST.C - Skeleton of program to make histogram of data from a log file. */

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


/* SKELETON OF PROGRAM TO MAKE HISTOGRAM FROM A LOG FILE.  This module is 
   linked with application-specific routines (the same as for plt.c) to 
   produce a program for building a histogram of values for a quantity
   obtained from a log file.  */

enum { PLT, TBL, HIST } program_type = HIST;


/* DISPLAY USAGE MESSAGE AND EXIT.  The plt_usage procedure must be defined 
   by the application specific module. */

void usage(void)
{
  extern void plt_usage(char *);

  plt_usage("hist quantities low high bins { log-file [ range ] }");

  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{ 
  char *quantities;
  double low, high, wid;
  int bins, frac;

  char **ap, **apa;

  int lindex, hindex, index_mod;

  log_file logf;
  log_gobbled logg;
  int very_first;

  quantities_described qd;
  quantities_held *qh;
  double *val;

  int total;
  int *bin;

  int N, V, a, i, b;

  /* Look at arguments other than log files and ranges. */

  if (argc<6) usage();

  quantities = argv[1];

  if ((low = atof(argv[2]))==0  && *argv[2]!='0'
   || (high = atof(argv[3]))==0 && *argv[3]!='0' || high<=low
   || (bins = atoi(argv[4]))==0)
  { usage();
  }

  frac = 0;
  if (bins<0)
  { frac = 1;
    bins = -bins;
  }

  wid = (high-low) / bins;

  ap = argv+5;

  for (apa = ap; *apa!=0 && strcmp(*apa,"/")!=0; apa++) ;

  if (ap==apa) usage();

  if (*apa!=0) apa += 1;

  for (a = 0; quant_app_arguments[a]; a++)
  { (*quant_app_arguments[a])(&apa);
  }

  if (*apa!=0) usage();

  /* Allocate space for bins, and set counts to zero. */

  bin = chk_alloc (bins, sizeof *bin);

  total = 0;

  for (i = 0; i<bins; i++) bin[i] = 0;

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

    N = quantities_requested(qd,quantities,1);
    quantities_available(qd,&logg);
    qh = quantities_storage(qd);

    V = 0;
    for (i = 0; i<N; i++) 
    { V += qd[i].high - qd[i].low + 1;
    }

    /* Go through all the records in the indicated range. */

    for (;;)
    {
      /* Get value for this iteration and increment appropriate bin. */
 
      quantities_evaluate(qd,qh,&logg);

      val = qh->value[0];

      for (i = 0; i<V; i++) 
      {
        if (val[i]>=low && val[i]<high)
        { b = (int) ( bins * ((val[i]-low)/(high-low)) );
          if (b<0 || b>=bins) abort();
          bin[b] += 1;
        }

        total += 1;
      }

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

    log_file_close(&logf);

    very_first = 0;
  }

  /* Output histogram. */

  for (b = 0; b<bins; b++)
  {
    if (frac)
    { printf( "%+.6e %.6f\n", low + wid*(b+0.5), (double)bin[b]/total);
    }
    else
    { printf( "%+.6e %6d\n", low + wid*(b+0.5), bin[b]);
    }
  }

  exit(0);
}
