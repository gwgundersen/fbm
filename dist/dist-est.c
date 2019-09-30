/* DIST-EST.C - Program to estimate the expectation of a function of state. */

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


/* VALUES AND IMPORTANCE WEIGHTS AT DATA POINTS. */

#define Max_points 200000		/* Maximum number of points allowed */

static int n_points;			/* Number of points used for estimate */

static double value[Max_points];	/* Value of formula at each point used*/
static float log_weight[Max_points];	/* Log of weight for each point used */

static double max_log_weight;		/* Maximum value of log weight */
static double w_mean;			/* Mean weight / exp(max_log_weight) */
static double ww_var;			/* Variance of normalized weights */
static double v_mean;			/* Weighted mean of formula value */
static double v_var;			/* Simple variance est. for value */

static double var_est_std;		/* Standard estimate for variance 
					   of v_mean */

static double v_mean_jack;		/* Jacknife estimate for mean of value*/
static double var_est_jack;		/* Jacknife estimate for variance 
                                           of v_mean_jack */

/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
    "Usage: dist-est formula [ temp-index ] { log-file [ range ] }\n");

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

  int first;

  char *value_formula;
  int tx;

  mc_temp_sched *sch;
  mc_temp_state *ts;
  mc_iter *it;

  dist_spec *dst;

  double *q;

  int nonzero_weight;

  double w, ww, jke, jkm;

  char **ap;

  char *a, *f;
  int c, i;

  /* Look at arguments other than log files and ranges. */

  if (argc<3) usage();

  value_formula = argv[1];

  if (strchr("0123456789",*argv[2])!=0)
  { tx = atoi(argv[2]);
    if (tx<0 || tx>=Max_temps) usage();
    ap = argv+3;
  }
  else
  { tx = -1;
    ap = argv+2;
  }

  /* Go through all the log files and ranges specified. */

  logg.req_size['i'] = sizeof (mc_iter);
  logg.req_size['o'] = sizeof (mc_ops);
  logg.req_size['t'] = sizeof (mc_traj);
  logg.req_size['b'] = sizeof (mc_temp_state);
  logg.req_size['m'] = sizeof (mc_temp_sched);
  logg.req_size['d'] = sizeof (dist_spec);

  n_points = 0;
  nonzero_weight = 0;
  first = 1;

  while (*ap!=0)
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

    /* Set things up the first time. */

    if (first)
    {
      if ((dst = logg.data['d'])==0)
      { fprintf(stderr,"No distribution specification in log file\n");
        exit(1);
      }

      a = dst->energy + strlen(dst->energy) + 1;
      if (dst->Bayesian) a += strlen(a) + 1;

      while (*a)
      { f = formula_def(a,&c,&i);
        formula_var[c][i] = formula(f,0,1,0);
        formula_var_exists[c][i] = 1;
        a += strlen(a) + 1;
      }

      (void) formula (dst->energy, 1, 0, 0);
      if (dst->Bayesian)
      { (void) formula (dst->energy + strlen(dst->energy) + 1, 1, 0, 0);
      }

      logg.req_size['q'] = dist_count_vars() * sizeof (double);
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
      sch = logg.data['m'];

      ts  = logg.data['b']!=0 && logg.index['b']==logg.last_index
             ? logg.data['b'] : 0;
      it  = logg.data['i']!=0 && logg.index['i']==logg.last_index
             ? logg.data['i'] : 0;
      q   = logg.data['q']!=0 && logg.index['q']==logg.last_index
             ? logg.data['q'] : 0;

      if (ts!=0 && sch==0)
      { fprintf(stderr,
          "Log file garbled: Has tempering state but no tempering schedule\n");
        exit(1);
      }

      /* Look at values here if tempering index is right. */

      if (ts==0 && tx<0 
       || ts!=0 && tx<0 && ts->inv_temp==1
       || ts!=0 && tx>=0 && ts->inv_temp==sch->sched[tx].inv_temp)
      { 
        if (q==0)
        { fprintf(stderr,"No variables stored with iteration\n");
          exit(1);
        }

        if (n_points>=Max_points)
        { fprintf(stderr,"Too many points for estimate (max %d)\n",Max_points);
          exit(1);
        }

        log_weight[n_points] = it ? it->log_weight : 0;
        if (log_weight[n_points]!=0) nonzero_weight = 1;

        if (n_points==0 || log_weight[n_points]>max_log_weight)
        { max_log_weight = log_weight[n_points];
        }

        dist_unpack_vars(q);
        value[n_points] = formula (value_formula, 0, 1, 0);

        n_points += 1;
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
 
  /* Find estimate, and other statistics. */

  if (n_points==0)
  { fprintf(stderr,"No records were found that could be used for estimate\n");
    exit(1);
  }

  w_mean = 0;

  for (i = 0; i<n_points; i++)
  { w = exp (log_weight[i] - max_log_weight);
    w_mean += w;
  }

  w_mean /= n_points;

  ww_var = 0;
  v_mean = 0;

  for (i = 0; i<n_points; i++)
  { ww = exp (log_weight[i] - max_log_weight) / w_mean;
    ww_var += (ww-1)*(ww-1);
    v_mean += ww * value[i];
  }

  ww_var /= n_points;
  v_mean /= n_points;

  var_est_std = 0;
  v_var = 0;

  for (i = 0; i<n_points; i++)
  { ww = exp (log_weight[i] - max_log_weight) / w_mean;
    v_var += ww * (value[i]-v_mean)*(value[i]-v_mean);
    var_est_std += ww*ww * (value[i]-v_mean)*(value[i]-v_mean);
  }

  v_var /= n_points;
  var_est_std /= (double)n_points*n_points;

  if (n_points>1)
  { 
    jkm = 0;

    for (i = 0; i<n_points; i++)
    { w = exp (log_weight[i] - max_log_weight);
      jke = (n_points*w_mean*v_mean - w*value[i]) / (n_points*w_mean - w);
      jkm += jke;
    }

    jkm /= n_points;

    v_mean_jack = n_points*v_mean - (n_points-1)*jkm;

    var_est_jack = 0;

    for (i = 0; i<n_points; i++)
    { w = exp (log_weight[i] - max_log_weight);
      jke = (n_points*w_mean*v_mean - w*value[i]) / (n_points*w_mean - w);
      var_est_jack += (jkm-jke)*(jkm-jke);
    }

    var_est_jack *= (double) (n_points-1) / n_points;
  }

  /* Print results. */

  printf("\nNumber of sample points: %d\n",n_points);
  if (nonzero_weight)
  { printf("\n  Variance of normalized weights: %g\n",ww_var);
    printf("  Adjusted sample size: %.1f  (reduction factor %.3f)\n",
      n_points/(1+ww_var), 1/(1+ww_var));
  }
  if (n_points/(1+ww_var)<10)
  { printf(
     "\nWARNING: Adjusted sample size < 10, results below may be unreliable\n");
  }

  if (nonzero_weight)
  { 
    printf("\nEstimates for importance weights:\n\n");
    printf("  Mean of weights: %g  (standard error %g)\n",
      w_mean * exp(max_log_weight), 
      w_mean * exp(max_log_weight) * sqrt(ww_var/n_points));
    printf("  Log of mean:    %+g  (standard error %g)\n",
      log(w_mean) + max_log_weight, sqrt(ww_var/n_points));
  }

  if (n_points>1 && nonzero_weight) 
  { printf("\nStandard estimates for %s:\n\n",value_formula);
  }
  else
  { printf("\nEstimates for %s:\n\n",value_formula);
  }
  printf("  Mean:    %g  (standard error %g)\n",v_mean,sqrt(var_est_std));
  printf("  Std.dev: %g\n",sqrt(v_var));
  if (nonzero_weight)
  { printf("\n  Effective sample size: %.1f  (reduction factor %.3f)\n",
      v_var/var_est_std, (v_var/var_est_std)/n_points);
  }

  if (n_points>1 && nonzero_weight)
  { printf("\nJacknife estimates for %s:\n\n",value_formula);
    printf("  Mean:    %g  (standard error %g)\n",
           v_mean_jack, sqrt(var_est_jack));
    if (nonzero_weight)
    { printf("\n  Effective sample size: %.1f  (reduction factor %.3f)\n",
        v_var/var_est_jack, (v_var/var_est_jack)/n_points);
    }
  }

  printf("\nNOTE: The standard errors assume points are independent\n\n");

  exit(0);
}
