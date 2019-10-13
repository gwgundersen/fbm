/* WMEAN.C - Program to estimate weighted means. */

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


/* VARIABLES RECORDING DATA. */

#define Max_batches 50000		/* Maximum number of batches allowed */
#define Max_values 20			/* Maximum number of values per line */

static int n_batches;			/* Number of batches */
static int n_values;			/* Number of values per line */

static double value[Max_batches][Max_values]; /* Mean values of each batch */
static double log_weight[Max_batches];	/* Log of weight for each batch */


/* STATISTICS COMPUTED. */

static double max_log_weight;		/* Maximum value of log batch weight */
static double w_mean;			/* Mean weight / exp(max_log_weight) */
static double ww_var;			/* Variance of normalized weights */
static double v_mean[Max_values];	/* Weighted means of values */
static double var_est[Max_values];	/* Estimates for variances of v_mean */


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr, "Usage: mean [ -W | -w ] [ -b ] [ batch-size ] <data\n");
  exit(1);
}


/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{ 
  int weighted, take_logs, bare, batch_size;
  double mx, w, v, vv[Max_values], ww;
  int i, n, k;
  char c;

  /* Look at arguments. */

  weighted = 0;
  take_logs = 0;
  bare = 0;
  batch_size = 1;

  for (;;)
  { if (argc>1 && strcmp(argv[1],"-w")==0)
    { weighted = 1;
      take_logs = 0;
    }
    else if (argc>1 && strcmp(argv[1],"-W")==0)
    { weighted = 1;
      take_logs = 1;
    }
    else if (argc>1 && strcmp(argv[1],"-b")==0)
    { bare = 1;
    }
    else
    { break;
    }
    argc -= 1;
    argv += 1;
  }

  if (argc>1)
  { if (sscanf(argv[1],"%d%c",&batch_size,&c)!=1 || batch_size<=0)
    { usage();
    }
  }

  if (argc>2) usage();

  /* Read weights and values from standard input, accumulating into 
     batch weights and batch means. */

  n_batches = 0;
  k = 0;

  for (n = 1; ; n++)
  { 
    scanf(" ");
    if (feof(stdin)) break;
    
    if (n_batches>=Max_batches)
    { fprintf(stderr,"Too many batches (max %d)\n",Max_batches);
      exit(1);
    }

    /* Read weight for line. */

    if (weighted)
    { if (scanf("%lf",&w)!=1 || take_logs && w<=0)
      { fprintf(stderr,"Bad weight on line %d\n",n);
        exit(1);
      }
    }
    else
    { w = 0;
    }

    /* Read values from line. */

    for (i = 0; ; i++)
    {
      while (scanf("%c",&c)==1 && c==' ') ;

      if (c=='\n') break;

      ungetc(c,stdin);
      if (scanf("%lf",&v)!=1)
      { fprintf(stderr,"Bad value on line %d\n",n);
        exit(1);
      }

      if (i>=Max_values)
      { fprintf(stderr,"Too many values per line (max %d)\n",Max_values);
        exit(1);
      }

      vv[i] = v;
    }

    /* Check that lines have the same number of values. */

    if (n_batches==0 && k==0) /* Very first line */
    { n_values = i;
    }
    else
    { if (i!=n_values)
      { fprintf(stderr,"Unequal numbers of values (line %d)\n",n);
        exit(1);
      }
    }

    /* Find log of weight if necessary.  Use very large negative log weight
       if weight is zero. */

    if (take_logs) 
    { w = w==0 ? -1e30 : log(w);
    }

    /* Combine lines making up a batch. */

    if (k==0) /* First line in batch */
    { log_weight[n_batches] = w;
      for (i = 0; i<n_values; i++) 
      { value[n_batches][i] = vv[i];
      }
    }
    else
    { if (weighted)
      { mx = log_weight[n_batches]>w ? log_weight[n_batches] : w;
        for (i = 0; i<n_values; i++)
        { value[n_batches][i] = 
           (exp(log_weight[n_batches]-mx)*value[n_batches][i] + exp(w-mx)*vv[i])
              / (exp(log_weight[n_batches]-mx) + exp(w-mx));
        }
        log_weight[n_batches] = mx+log(exp(log_weight[n_batches]-mx)+exp(w-mx));
      }
      else
      { for (i = 0; i<n_values; i++)
        { value[n_batches][i] = (k*value[n_batches][i] + vv[i]) / (k+1);
        }
        log_weight[n_batches] = log((double)(k+1));
      }
    }

    k += 1;

    /* Handle end of batch. */

    if (k==batch_size)
    {
      if (n_batches==0 || log_weight[n_batches]>max_log_weight)
      { max_log_weight = log_weight[n_batches];
      }

      n_batches += 1;
      k = 0;
    }
  }

  if (k!=0)
  { fprintf (stderr,
      "WARNING:  Last batch is incomplete (%d out of %d) - ignored\n",
      k, batch_size);
  }

  if (n_batches==0)
  { fprintf(stderr,"No batches with non-zero weight\n");
    exit(1);
  }

  if (0) /* For debugging */
  { for (n = 0; n<n_batches; n++)
    { printf("%.3f",log_weight[n]);
      for (i = 0; i<n_values; i++)
      { printf(" %.3f",value[n][i]);
      }
      printf("\n");
    }
  }

  /* Compute statistics on weights. */

  if (weighted)
  { w_mean = 0;
    for (n = 0; n<n_batches; n++)
    { w = exp (log_weight[n] - max_log_weight);
      w_mean += w;
    }
    w_mean /= n_batches;
  }
  else
  { w_mean = 1;
  }

  ww_var = 0;
    
  if (weighted && n_batches>1)
  { for (n = 0; n<n_batches; n++)
    { ww = exp (log_weight[n] - max_log_weight) / w_mean;
      ww_var += (ww-1)*(ww-1);
    }
    ww_var /= n_batches-1;
  }

  /* Compute weighted means for each value. */

  for (i = 0; i<n_values; i++)
  { v_mean[i] = 0;
  }

  for (n = 0; n<n_batches; n++)
  { if (weighted)
    { ww = exp (log_weight[n] - max_log_weight) / w_mean;
      for (i = 0; i<n_values; i++)
      { v_mean[i] += ww * value[n][i];
      }
    }
    else
    { for (i = 0; i<n_values; i++)
      { v_mean[i] += value[n][i];
      }
    }
  }

  for (i = 0; i<n_values; i++)
  { v_mean[i] /= n_batches;
  }

  /* Find standard errors for means of values. */

  for (i = 0; i<n_values; i++)
  { var_est[i] = 0;
  }

  for (n = 0; n<n_batches; n++)
  { ww = exp (log_weight[n] - max_log_weight) / w_mean;
    for (i = 0; i<n_values; i++)
    { var_est[i] += ww*ww * (value[n][i]-v_mean[i])*(value[n][i]-v_mean[i]);
    }
  }

  for (i = 0; i<n_values; i++)
  { var_est[i] /= (double)n_batches*n_batches;
  }

  /* Print results. */

  if (!bare)
  { printf("\nNumber of batches: %d\n",n_batches);
    printf("Number of lines per batch: %d\n",batch_size);

    if (weighted && n_batches>1)
    { printf("\nVariance of normalized weights: %g\n",ww_var);
      printf("Adjusted sample size: %.1f  (reduction factor %.3f)\n",
        n_batches/(1+ww_var), 1/(1+ww_var));
      if (n_batches/(1+ww_var)<10)
      { printf(
     "\nWARNING: Adjusted sample size < 10, results below may be unreliable\n");
      }
    }
  }

  if (weighted)
  { if (bare)
    { if (n_values==0)
      { printf ("%+.13e %.13e\n",
         log(w_mean/batch_size) + max_log_weight, sqrt(ww_var/n_batches));
      }
    }
    else
    { printf("\nLog of mean weight: %+.2f (standard error %.2f)\n",
       log(w_mean/batch_size) + max_log_weight, sqrt(ww_var/n_batches));
    }
  }

  if (!bare && n_values>0)
  { printf ("\n         Mean        Standard error\n\n");
  }

  for (i = 0; i<n_values; i++)
  { if (bare)
    { printf("%+.13e %.13e\n",v_mean[i],sqrt(var_est[i]));
    }
    else
    { printf("%16.7f  %14.7f\n",v_mean[i],sqrt(var_est[i]));
    }
  }

  if (!bare) 
  { printf("\n");
  }

  exit(0);
}
