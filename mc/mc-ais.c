/* MC-AIS.C - Program to show how well Annealed Importance Sampling worked. */

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


/* DATA FOR EACH TEMPERATURE. */

static double w_mean[Max_temps];	/* Mean of weights (as log) */
static double ww_var[Max_temps];	/* Variance of normalized weights */
static double logw_mean[Max_temps];	/* Mean of log of weights */
static double logw_var[Max_temps];	/* Variance of log of weights */
static int    n_points[Max_temps];	/* Number of points */


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,"Usage: mc-ais values { log-file [ range ] }\n");
  fprintf(stderr,
   "Values: i=inv-temp T=temp I=temp-index m=mean-weight M=log(m) F=-log(m)\n");
  fprintf(stderr,
   "  v=var(norm-wts) V=var(log-wts) W=log(1+var(norm-wts)) a=adj-sample-size\n");
  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{ 
  int lindex, hindex, index_mod;

  log_file logf;
  log_gobbled logg;

  char *values;

  int first, pass;

  mc_temp_sched *sch;
  mc_temp_state *ts;
  mc_iter *it;

  char **ap;

  double ww;
  char *v;
  int i;

  /* Look at arguments other than log files and ranges. */

  if (argc<3) usage();

  values = argv[1];

  /* Set records sizes required. */

  logg.req_size['i'] = sizeof (mc_iter);
  logg.req_size['o'] = sizeof (mc_ops);
  logg.req_size['t'] = sizeof (mc_traj);
  logg.req_size['b'] = sizeof (mc_temp_state);
  logg.req_size['m'] = sizeof (mc_temp_sched);

  /* Go twice through all the log files and ranges specified.  Compute
     means the first time, variances the second. */

  for (pass = 1; pass<=2; pass++)
  {
    first = 1;

    ap = argv+2;  

    while (*ap)
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
  
      log_gobble_init(&logg,!first);
      
      /* Gobble up records with negative indexes. */
  
      while (logf.header.index<0)
      { log_gobble(&logf,&logg);
      }
  
      /* Skip to start of range, gobble up records at start. */
  
      while (!logf.at_end 
         && (logf.header.index<lindex || logf.header.index%index_mod!=0))
      { log_file_forward(&logf);
      }
  
      if (logf.at_end) continue;
  
      log_gobble(&logf,&logg);
  
      /* Go through all the records in the indicated range. */
  
      for (;;)
      {
        /* See what's here. */
  
        sch = logg.data['m'];
  
        if (sch==0)
        { fprintf(stderr,"No tempering schedule present\n");
          exit(1);
        }
  
        ts  = logg.data['b']!=0 && logg.index['b']==logg.last_index
               ? logg.data['b'] : 0;
        it  = logg.data['i']!=0 && logg.index['i']==logg.last_index
               ? logg.data['i'] : 0;
  
        /* Look at data for this index, if anything is happening. */
  
        if (ts!=0 || it!=0)
        {
          if (ts==0)
          { fprintf(stderr,"No tempering state present (%s %d)\n",
               logf.file_name, logg.last_index);
            exit(1);
          }
  
          if (it==0)
          { fprintf(stderr,"No record describing iteration present (%s %d)\n",
               logf.file_name, logg.last_index);
            exit(1);
          }
  
          i = mc_temp_index(sch,ts->inv_temp);
  
          /* Look at data for pass 1. */

          if (pass==1)
          {
            if (n_points[i]==0) 
            { w_mean[i] = it->log_weight;
              logw_mean[i] = it->log_weight;
            }
            else
            { w_mean[i] = addlogs (w_mean[i], it->log_weight);
              logw_mean[i] = logw_mean[i] + it->log_weight;
            }
  
            n_points[i] += 1;
          }

          /* Look at data for pass 2. */

          if (pass==2)
          { 
            ww = exp(it->log_weight-w_mean[i]);
            ww_var[i] += (ww-1)*(ww-1);
            logw_var[i] += (it->log_weight - logw_mean[i])
                             * (it->log_weight - logw_mean[i]);
          }
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
      first = 0;
    }

    /* Finish up pass 1. */

    if (pass==1)
    { for (i = 0; i<Max_temps; i++)
      { if (n_points[i]>0)
        { w_mean[i] -= log((double)n_points[i]);
          logw_mean[i] /= n_points[i];
        }
        if (sch->sched[i].inv_temp==1) break;
      }
    }

    /* Finish up pass 2. */

    if (pass==2)
    { for (i = 0; i<Max_temps; i++)
      { if (n_points[i]>1)
        { ww_var[i] /= (n_points[i]-1);
          logw_var[i] /= (n_points[i]-1);
        }
        if (sch->sched[i].inv_temp==1) break;
      }
    }
  }
  
  /* Print values for each temperature. */

  for (i = 0; i<Max_temps; i++)
  {
    if (n_points[i]>0)
    { 
      for (v = values; *v; v++)
      { 
        switch (*v)
        { 
          case 'I': 
          { printf(" %4d",i);
            break;
          }
  
          case 'i':
          { printf(" %.6f",sch->sched[i].inv_temp);
            break;
          }
   
          case 'T':
          { printf(" %9.4f",1/sch->sched[i].inv_temp);
            break;
          }
  
          case 'm':
          { printf(" %.6e",
              n_points[i]==0 ? 1 : exp(w_mean[i]));
            break;
          }
  
          case 'M':
          { printf(" %.6e",
              n_points[i]==0 ? 0 : w_mean[i]);
            break;
          }
  
          case 'F':
          { printf(" %.6e",
              n_points[i]==0 ? 0 : -w_mean[i]);
            break;
          }

          case 'v':
          { printf(" %8.3f", n_points[i]>1 ? ww_var[i] : 0);
            break;
          }

          case 'V':
          { printf(" %8.3f", n_points[i]>1 ? logw_var[i] : 0);
            break;
          }

          case 'a':
          { printf(" %9.3f", n_points[i]>1 ? n_points[i]/(1+ww_var[i]) : 0);
            break;
          }

          case 'W':
          { printf(" %9.3f", n_points[i]>1 ? log (1+ww_var[i]) : 0);
            break;
          }
  
          default:
          { fprintf(stderr,"Invalid value specifier: %c\n",*v);
            exit(1);
          }       
        }
      }
  
      printf("\n");
    }

    if (sch->sched[i].inv_temp==1) break;
  }
 
  exit(0);
}


/* FIND INDEX FOR SIMULATED TEMPERING INVERSE TEMPERATURE.  Copied from
   mc-util.c, since unfortunately other routines in that file reference
   application-specific routines. */

int mc_temp_index
( mc_temp_sched *sch,		/* Schedule of inverse temperatures */
  float inv_temp		/* Inverse temperature to find */
)
{ 
  int i;

  if (sch==0)
  { fprintf(stderr,"No tempering schedule has been specified\n");
    exit(1);
  }

  for (i = 0; sch->sched[i].inv_temp!=inv_temp; i++)
  { if (i==Max_temps-1) abort();
  }

  return i;
}
