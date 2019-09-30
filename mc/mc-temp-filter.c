/* MC-TEMP-FILTER.C - Program to filter out iterations at a given temperature.*/

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


static void usage(void);


/* FIND INDEX FOR SIMULATED TEMPERING INVERSE TEMPERATURE.  Extracted
   from mc-util because of linking problems. */

static int temp_index
( mc_temp_sched *sch,		/* Schedule of inverse temperatures */
  float inv_temp		/* Inverse temperature to find */
)
{ 
  int i;

  if (sch==0) abort();

  for (i = 0; sch->sched[i].inv_temp!=inv_temp; i++)
  { if (i==Max_temps-1) abort();
  }

  return i;
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  log_file logf_in, logf_out;
  log_gobbled logg;

  int seq, found_highest, tours;
  int index;

  mc_temp_sched *sch;
  mc_temp_state *ts;
  float *temp;
  int i;

  if (argc!=3 && argc!=4) usage();

  index = 0;

  if (argc==4)
  { index = atoi(argv[1]);
    if (index<=-Max_temps || index>=Max_temps) 
    { fprintf(stderr,"Temperature index out of range\n");
      exit(1);
    }
  }

  logf_in.file_name = argv[argc-2];
  log_file_open(&logf_in,0);

  logf_out.file_name = argv[argc-1];
  log_file_create(&logf_out);

  log_gobble_init(&logg,0);

  logg.req_size['b'] = sizeof (mc_temp_state);
  logg.req_size['m'] = sizeof (mc_temp_sched);

  sch = 0;
  tours = 0;
  found_highest = 0;
  seq = 1;

  while (!logf_in.at_end)
  {
    log_gobble(&logf_in,&logg);

    ts = logg.data['b'];
    sch = logg.data['m'];

    if (argc==3 && index==0) index = temp_index(sch,1.0);

    temp = ts==0 || logg.index['b']!=logg.last_index ? 0 : &ts->inv_temp;

    if (found_highest && temp!=0 && *temp==1)
    { tours += 1;
      found_highest = 0;
    }

    if (sch!=0 && temp!=0 && *temp==sch->sched[0].inv_temp)
    { found_highest = 1;
    }

    if (logg.last_index<0 
     || temp!=0 && index<0 
        && (*temp==sch->sched[0].inv_temp || *temp==sch->sched[-index].inv_temp)
     || temp!=0 && index>=0 && *temp==sch->sched[index].inv_temp)
    {
      for (i = 0; i<128; i++)
      { 
        if (logg.data[i]!=0 && logg.index[i]==logg.last_index)
        { 
          logf_out.header.index = logg.last_index>0 ? seq : logg.index[i];
          logf_out.header.type  = i;
          logf_out.header.size  = logg.actual_size[i];
          log_file_append (&logf_out, logg.data[i]);
        }
      }

      if (logg.last_index>0) seq += 1;
    }

  }

  log_file_close(&logf_in);
  log_file_close(&logf_out);

  if (sch!=0)
  { fprintf(stderr,"Number of tours to highest temperature and back: %d\n",
                    tours);
  }
  
  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */

void usage(void)
{ 
  fprintf(stderr, 
   "Usage: mc-temp-filter [-]temp-index input-log-file output-log-file\n");
  exit(1);
}
