/* TBL.C - Skeleton of program to output table of quantities from log files. */

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


/* SKELETON OF PROGRAM TO OUTPUT TABLE OF DATA FROM LOG FILES.  This module is 
   linked with application-specific routines to produce a program for writing
   quantities obtained from a log file to standard output, in a tabular form.
   This form is suitable for viewing, or for reading into S-Plus. */

enum { PLT, TBL, HIST } program_type = TBL;


/* DISPLAY USAGE MESSAGE AND EXIT.  The plt_usage procedure must be defined 
   by the application specific module. */

void usage(void)
{
  extern void plt_usage(char *);

  plt_usage("tbl [ -h ] quantities { log-file [ range ] }");

  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{ 
  char *quant;
  char **ap, **apa;
  int header;

  int lindex, hindex, index_mod;

  log_file logf;
  log_gobbled logg;

  int very_first;

  quantities_described qd, qd_first;
  quantities_held *qh_first, *qh, **qhp;

  char str[20];
  int a, i, v;
  int N, V;

  /* Look at flag arguments. */

  header = 0;
 
  while (argc>1 && argv[1][0]=='-')
  {
    if (strcmp(argv[1],"-h")==0) 
    { header = 1;
    }
    else 
    { usage();
    }

    argc -= 1;
    argv += 1;
  }

  /* Look at other arguments before log files and ranges. */

  if (argc<3) usage();

  quant = argv[1];

  ap = argv+2;

  /* Look at application-specific arguments. */

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

    N = quantities_requested(qd,quant,1);

    quantities_available(qd,&logg);

    if (very_first)
    { for (i = 0; i<Max_quantities; i++)
      { qd_first[i] = qd[i];
      }
    }
    else
    { for (i = 0; i<Max_quantities; i++)
      { if (qd[i].low!=qd_first[i].low || qd[i].high!=qd_first[i].high)
        { fprintf(stderr,
           "Range of values for '%c' quantity differs between log files\n",
            qd[i].letter);
          exit(1);
        }
      }
    }

    V = 0;

    for (i = 0; i<N; i++)
    { V += qd[i].high - qd[i].low + 1;
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

    /* Print headers, if desired and this is the first file. */

    if (very_first && header)
    { for (i = 0; i<N; i++)
      { for (v = 0; v<qd[i].high-qd[i].low+1; v++)
        { sprintf(str,"%c",qd[i].letter);
          if (qd[i].modifier!=-1) 
          { sprintf(str+strlen(str),"%d",qd[i].modifier);
          }
          if (qd[i].low!=-1)
          { sprintf(str+strlen(str),"@%d",qd[i].low+v);
          }
          printf("%11s     ",str);
        }
      }
      printf("\n");
    }

    /* Plot data collected. */

    for (qh = qh_first; qh!=0; qh = qh->next)
    {
      for (v = 0; v<V; v++)
      { printf ("%15.8e ", *(qh->value[0]+v));
      }

      printf("\n");
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

    very_first = 0;
  }

  exit(0);
}
