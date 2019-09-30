/* MC-TEMP-SCHED.C - Specify tempering schedule for Markov chain simulation.*/

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


static void usage (void);
static void display_sched (log_gobbled *);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  static mc_temp_sched sch0, *sch = &sch0; 

  float inv_temp, bias;
  int no_sched;

  log_file logf;
  log_gobbled logg;

  int arithmetic;
  char **ap;
  char *p;
  int t, n, i;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing schedule. */

  if (argc==2)
  {
    no_sched = 1;

    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    logg.req_size['m'] = sizeof *sch;

    if (!logf.at_end) 
    { do 
      { log_gobble(&logf,&logg); 
      } while (!logf.at_end && logg.last_index<-1);
  
      if (logg.data['m']!=0)
      { printf("\nTEMPERING SCHEDULE IN LOG FILE:\n");
        display_sched(&logg);
        no_sched = 0;
      }
    }

    if (no_sched)
    { printf("\nNo tempering schedule found\n\n");
    }

    exit(0);
  }

  /* Handle null schedule. */

  if (argc==3 && strcmp(argv[2],"-")==0)
  { argc -= 1;
  }

  /* Otherwise, start by examining temperature/bias arguments in reverse
     order, storing them the wrong way around in the schedule. */

  sch->sched[0].inv_temp = 1;
  sch->sched[0].bias = 0;
  t = 1;

  for (ap = argv+argc-1; ap>argv+1; ap--)
  { 
    arithmetic = 0;
    n = 1;

    if (strchr(*ap,':')!=0)
    { p = strchr(*ap,':')+1;
      if (*p=='+') 
      { arithmetic = 1;
        p += 1;
      }
      if (*p==0) usage();
      n = atoi(p);
      if (n<=0) usage();
    }

    if (t+n>Max_temps) 
    { fprintf(stderr,"Too many temperatures in schedule (max %d)\n",
              Max_temps-1);
      exit(1);
    }

    if (**ap=='/')
    {
      inv_temp = atof(*ap+1);
      if (inv_temp<=1)
      { fprintf(stderr,"Temperatures must be greater than one\n");
        exit(1);
      }
      inv_temp = 1 / inv_temp;
    }
    else
    {
      inv_temp = atof(*ap);
      if (inv_temp<0 || inv_temp>=1)
      { fprintf(stderr,"Inverse temperatures must be in the range [0,1)\n");
        exit(1);
      }
    }

    if (strchr(*ap,'@')!=0)
    { bias = atof(strchr(*ap,'@')+1);
    }
    else
    { bias = 0;
    }

    if (inv_temp==0 && n>1 && !arithmetic)
    { fprintf(stderr,
        "Can't start geometric series with inverse temperature of 0\n");
      exit(1);
    }

    if (n==1)
    { sch->sched[t].inv_temp = inv_temp;
      sch->sched[t].bias = bias;
    }
    else
    { for (i = n-1; i>=0; i--)
      { if (arithmetic)
        { sch->sched[t+n-1-i].inv_temp = 
            inv_temp + (sch->sched[t-1].inv_temp-inv_temp)*i/n;
        }
        else
        { sch->sched[t+n-1-i].inv_temp = exp 
            (log(inv_temp)+(log(sch->sched[t-1].inv_temp)-log(inv_temp))*i/n);
        }
        sch->sched[t+n-1-i].bias = 
          bias + (sch->sched[t-1].bias-bias)*i/n;
      }
    }

    t += n;
  }

  /* Reverse the schedule. */

  for (i = 0; i<t/2; i++)
  { struct mc_temp_bias x;
    x = sch->sched[i];
    sch->sched[i] = sch->sched[t-i-1];
    sch->sched[t-i-1] = x;
  }

  /* Write schedule to log file. */

  log_file_open (&logf, 1);

  logf.header.type = 'm';
  logf.header.index = -1;
  logf.header.size = sizeof *sch;
  log_file_append (&logf, sch);

  log_file_close (&logf);

  exit(0);
}


/* DISPLAY SCHEDULE. */

static void display_sched
( log_gobbled *logg
)
{
  mc_temp_sched *sch;
  int t;

  sch = logg->data['m'];

  printf("\nIndex  Inv-temp     Bias  Temp\n\n");

  for (t = 0; t<Max_temps; t++)
  { printf ("%4d   %.6f  @%+6.1f", 
      t, sch->sched[t].inv_temp, sch->sched[t].bias);
    if (sch->sched[t].inv_temp==0)
    { printf("  /inf\n");
    }
    else
    { printf("  /%.1f\n", 1.0/sch->sched[t].inv_temp);
    }
    if (sch->sched[t].inv_temp==1) break;
  }

  printf("\n");
  
  if (t==Max_temps)
  { printf("WARNING: No temperature of one at end of schedule!\n");
  }
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr, 
    "Usage: mc-temp_sched log-file { inv-temperature[@bias][:[+]n] }\n");

  fprintf(stderr, 
    "   or: mc-temp_sched log-file -  (to specify a null schedule)\n");

  fprintf(stderr,
    "   or: mc-temp-sched log-file    (to display stored schedule)\n");

  exit(1);
}
