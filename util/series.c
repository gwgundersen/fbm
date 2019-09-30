/* SERIES.C - Program to analyse stationary time series data. */

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


#define Max_realizations   200
#define Max_length       50000

static double *data[Max_realizations];	/* Data for each realization */
static int length[Max_realizations];	/* Length of each realization */
static int nr;				/* Number of realizations */

static double cvs[Max_length];		/* Sum for calculating covariances */
static double cvn[Max_length];		/* Number of items in each sum */
static int ml;				/* Maximum lag with any data */

static void usage(void);


main
( int argc,
  char **argv
)
{ 
  int op_m, op_s, op_v, op_a, op_c, op_b, op_e;
  char *options;
  int max_lag, have_presumed_mean;
  double presumed_mean;

  double mean, submean[Max_realizations], subvar[Max_realizations];
  int same_lengths;
  double variance;

  int c, r, t, l, n;
  double d, s;

  /* Look at arguments. */

  if (argc<2 || argc>4) usage();

  options = argv[1];

  op_m = strchr(options,'m')!=0;
  op_s = strchr(options,'s')!=0;
  op_v = strchr(options,'v')!=0;
  op_a = strchr(options,'a')!=0;
  op_c = strchr(options,'c')!=0;
  op_b = strchr(options,'b')!=0;
  op_e = strchr(options,'e')!=0;

  if (strlen(options) != op_m+op_s+op_v+op_a+op_c+op_b+op_e) usage();

  if (op_b && !op_a && !op_c)
  { fprintf(stderr,"Option b makes no sense without a or c\n");
    exit(1);
  }

  max_lag = -1;
  have_presumed_mean = 0;

  if (argc>2)
  { max_lag = atoi(argv[2]);
    if (max_lag<=0 && strcmp(argv[2],"0")!=0) usage();
  }

  if (argc>3)
  { have_presumed_mean = 1;
    presumed_mean = atof(argv[3]);
  }

  if (op_m && have_presumed_mean)
  { fprintf(stderr,"Option m makes no sense with a presumed mean\n");
    exit(1);
  }

  if ((op_m || op_v || op_s) && op_b)
  { fprintf(stderr,"Options m, v, and s make no sense with option b\n");
    exit(1);
  }

  /* Read data. */

  nr = 0;
 
  for (;;)
  { 
    /* Skip to start of next realization. */

    do { c = getc(stdin); } while (c==' ' || c=='\t' || c=='\n');

    if (c==EOF) break;

    if (nr==Max_realizations) 
    { fprintf(stderr,"Too many realizations of series (max %d)\n",
              Max_realizations);
      exit(1);
    }

    /* Read data for this realization. */

    data[nr] = chk_alloc (Max_length, sizeof (double));

    ungetc (c, stdin);

    length[nr] = 0;

    for (;;)
    { 
      if (scanf("%lf",&d)!=1)
      { fprintf(stderr,"Bad data in time series\n");
        exit(1);
      }

      if (length[nr]==Max_length)
      { fprintf(stderr,"Realization of time series is too long (max %d)\n",
                Max_length);
        exit(1);
      }

      data[nr][length[nr]] = d;
      length[nr] += 1;

      do { c = getc(stdin); } while (c!=EOF && c!='\n');

      if (c==EOF) break;
   
      do { c = getc(stdin); } while (c==' ' || c=='\t');

      if (c=='\n' || c==EOF) break;

      ungetc (c, stdin);
    }

    nr += 1;

    if (c==EOF) break;
  }

  if (nr==0) exit(1);

  /* Just echo data, if that's what we're asked to do. */

  if (op_e)
  {
    for (r = 0; r<nr; r++)
    {
      if (r>0) printf("\n");
 
      for (t = 0; t<length[r]; t++)
      { printf("%20.8e\n",data[r][t]);
      }
    }

    exit(0);
  }

  /* Calculate mean and submeans. */

  s = 0;
  n = 0;

  same_lengths = 1;

  for (r = 0; r<nr; r++)
  { submean[r] = 0;
    for (t = 0; t<length[r]; t++)
    { s += data[r][t];
      submean[r] += data[r][t];
      n += 1;
    }
    submean[r] /= length[r];
    if (r>0 && length[r-1]!=length[r]) same_lengths = 0;
  }

  mean = have_presumed_mean ? presumed_mean : s/n;

  /* Calculate variance. */

  s = 0;
  n = 0;

  for (r = 0; r<nr; r++)
  { for (t = 0; t<length[r]; t++)
    { d = data[r][t] - mean;
      s += d*d;
      n += 1;
    }
  }

  variance = s/n;

  /* Calculate autocovariances. */

  ml = -1;
  
  for (r = 0; r<nr; r++)
  {
    for (t = 0; t<length[r]; t++)
    { for (l = 0; t+l<length[r]; l++)
      { if (max_lag>=0 && l>max_lag) break;
        if (l>ml)
        { cvs[l] = 0;
          cvn[l] = 0;
          ml = l;
        }
        cvs[l] += (data[r][t] - mean) * (data[r][t+l] - mean);
        cvn[l] += 1;
      }
    }
  }

  /* Print counts, if not giving bare listings. */

  if (!op_b)
  { 
    n = 0;
    for (r = 0; r<nr; r++) n += length[r];

    printf("\nNumber of realizations: %d  Total points: %d\n",nr,n);
  }

  /* Print presumed mean, if there is one, and bare listing not desired. */

  if (!op_b && have_presumed_mean)
  { printf("\nPresumed mean: %.6f\n", presumed_mean);
  }

  /* Print simple statistics, if desired. */

  if (op_m)
  {
    printf("\n");

    printf("Mean: %.6f  ",mean);

    /* Print standard error estimated from multiple realizations. */

    if (nr>1 && same_lengths)
    { s = 0;
      for (r = 0; r<nr; r++) 
      { d = submean[r] - mean;
        s += d*d;
      }
      s /= nr-1;
      printf("S.E. from submeans: %.6f  ",sqrt(s/nr));
    }

    /* Print standard error estimated from correlations. */

    if (max_lag>=0)
    { 
      for (r = 0; r<nr; r++)
      { subvar[r] = variance;
        for (l = 1; l<=ml && l<length[r]; l++)
        { subvar[r] += 2 * (1 - (double)l/length[r]) * (cvs[l]/cvn[l]);
        }
        subvar[r] /= length[r];
      }
 
      s = 0;

      for (r = 0; r<nr; r++) 
      { s += subvar[r] * (length[r]*length[r]) / (n*n);
      }

      printf("S.E. from correlations: %.6f",sqrt(s));
    }

    printf("\n");
  }

  if (op_v || op_s)
  { printf("\n");
    if (op_v) printf("Variance: %.6f  ",variance);
    if (op_s) printf("Standard deviation: %.6f  ",sqrt(variance));
    printf("\n");
  }

  /* Print autocorrelation data, if desired. */

  if ((op_a || op_c) && variance>0)
  {
    if (!op_b)
    { printf("\n  Lag  ");
      if (op_a) printf("Autocorr.  ");
      if (op_c) printf("Cum. Corr.");
      printf("\n\n");
    }

    s = 1;

    for (l = 1; l<=ml; l++)
    {
      printf("%5d  ",l);
      
      if (op_a) printf("%9.6f ",(cvs[l]/cvn[l]) / variance);

      s += 2 * (cvs[l]/cvn[l]) / variance;

      if (op_c) printf("%11.6f",s);

      printf("\n");
    }
  }

  exit(0);
}

static void usage(void)
{ 
  fprintf(stderr,"Usage: series options [ max-lag [ presumed-mean ] ] <data\n");
  exit(1);
}
