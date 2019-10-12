/* SRC-INTENSITY.C - Make predictions for source intensity in grid cells.

/* Copyright (c) 2007 by Radford M. Neal 
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
#include "rand.h"
#include "log.h"
#include "mc.h"
#include "src.h"


#define DEBUG 0			/* Normally 0 */


/* LOCAL PROCEDURES. */

static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  double w, f, lf, sum_weights, max_log_weight;
  int lindex, hindex, mod_no;
  double target_index;
  int nx, ny, nz, nt;
  int hx, hy, hz, ht;
  log_gobbled logg;
  log_file logfile;
  int N_records_used;
  mc_temp_state *ts;
  mc_iter *it;
  char junk;
  char **ap;

  src_spec *src;
  src_params *params;
  double highest_time;
  int xi, yi, zi, ti;
  float *cells;
  double Q;
  int i, j;

  /* Find number of cells in each dimension. */

  nx = ny = nz = nt = 1;
  hx = hy = hz = ht = 0;

  if (argc<2) usage();

  ap = argv+1;

  if (*ap!=0 && strcmp(*ap,"/")!=0)
  { if (sscanf(*ap,"%d%c",&nx,&junk)!=1 || nx<1) usage();
    hx = 1;
    ap += 1;
  }

  if (*ap!=0 && strcmp(*ap,"/")!=0)
  { if (sscanf(*ap,"%d%c",&ny,&junk)!=1 || nx<1) usage();
    hy = 1;
    ap += 1;
  }

  if (*ap!=0 && strcmp(*ap,"/")!=0)
  { if (sscanf(*ap,"%d%c",&nz,&junk)!=1 || nx<1) usage();
    hz = 1;
    ap += 1;
  }

  if (*ap!=0 && strcmp(*ap,"/")!=0)
  { if (sscanf(*ap,"%d%c",&nt,&junk)!=1 || nx<1) usage();
    ht = 1;
    ap += 1;
  }

  if (strcmp(*ap,"/")!=0) usage();
  ap += 1;

  if (DEBUG) fprintf(stderr,"nx %d ny %d nz %d nt %d\n",nx,ny,nz,nt);

  /* Allocate space for cells. */

  cells = chk_alloc (nx*ny*nz*nt, sizeof *cells);
  for (i = nx*ny*nz*nt-1; i>=0; i--)
  { cells [i] = 0;
  }

  /* Open first log file and get specifications. */

  if (*ap==0) usage();

  logfile.file_name = *ap++;

  log_file_open (&logfile, 0);

  log_gobble_init(&logg,0);
  src_record_sizes(&logg);

  while (!logfile.at_end && logfile.header.index<0)
  { log_gobble(&logfile,&logg);
  }

  src = (src_spec *) logg.data['S'];
  if (src==0)
  { fprintf(stderr,"No source specification in log file\n");
    exit(1);
  }

  if (nx>1 && src->low[0]==src->high[0]
   || nx>1 && src->low[0]==src->high[0]
   || nx>1 && src->low[0]==src->high[0]
   || nt>1 && src->max_stop==0)
  { fprintf (stderr,
       "More than one cell is not allowed for dimensions with zero size\n");
    exit(1);
  }

  highest_time = src->max_start + src->max_duration;
  if (src->max_stop<highest_time)
  { highest_time = src->max_stop;
  }

  if (nt>1 && highest_time>=1e30)
  { fprintf (stderr, 
       "More than one time cell is not allowed when sources don't stop\n");
    exit(1);
  }

  /* Go through all the iterations, taken from all the log files. */

  rand_seed(1);

  N_records_used = 0;
  sum_weights =  0;

  for (;;)
  {
    /* Figure out range to use from this log file. */

    if (*ap==0) usage();

    parse_range (*ap++, &lindex, &hindex, &mod_no);

    if (hindex==-2) 
    { hindex = lindex;
    }
    else if (hindex==-1) 
    { log_file_last(&logfile);
      hindex = logfile.header.index;
      log_file_first(&logfile); 
    }

    /* Reduce number of iterations asked for if it's more than is available.*/

    if (mod_no<0 && hindex-lindex+1<-mod_no)
    { mod_no = - (hindex-lindex+1);
    }

    if (mod_no<0) target_index = lindex;

    for (;;)
    {
      /* Skip to next desired index, or to end of range. */

      while (!logfile.at_end && logfile.header.index<=hindex
       && (logfile.header.index<lindex || 
           (mod_no>0 ? logfile.header.index%mod_no!=0
                     : logfile.header.index<(int)(target_index+0.5))))
      { 
        log_file_forward(&logfile);
      }

      if (logfile.at_end || logfile.header.index>hindex)
      { break;
      }
  
      /* Gobble up records for this index. */

      log_gobble(&logfile,&logg);

      if (DEBUG)
      { fprintf (stderr, "Using record from log file %s at index %d.\n",
          logfile.file_name, logg.last_index);
      }

      /* Look at weight in iteration record and tempering state. */

      it = logg.data['i'] != 0 && logg.index['i']==logg.last_index
             ? logg.data['i'] : 0;
      ts  = logg.data['b']!=0 && logg.index['b']==logg.last_index
             ? logg.data['b'] : 0;

      w = it!=0 && it->log_weight!=0 ? it->log_weight : 0;
     
      /* Use records at this index to make predictions for the test cases. */

      if (ts!=0 && ts->inv_temp!=1)
      { /* Ignore - wrong temperature */
      }
      else if (logg.index['q']!=logg.last_index)
      { fprintf(stderr,
          "Warning: Missing data at index %d in %s - ignored this index\n",
           logg.last_index, logfile.file_name);
      }
      else
      { 
        /* Adjust weights as necessary, accounting for new maximum. */

        if (N_records_used==0)
        { max_log_weight = w;
          lf = 0;
          f = 1;
        }
        else if (w>max_log_weight)
        { lf = max_log_weight-w;
          f = exp(lf);
          sum_weights *= f;
          for (i = nx*ny*nz*nt-1; i>=0; i--)
          { cells[i] *= f;
          }
          max_log_weight = w;
          lf = 0;
          f = 1;
        }
        else
        { lf = w-max_log_weight;
          f = exp(lf);
        }

        /* Add weighted contribution to cells where the sources in this
           iteration are located. */

        params = (src_params *) logg.data['q'];

        sum_weights += f;

        for (i = 0; i<(int)params->N0; i++)
        { 
          xi = (int) (nx * (params->src[i].coord[0] - src->low[0]) 
                         / (src->high[0] - src->low[0]));
          if (xi==nx) xi = nx-1;
          yi = (int) (ny * (params->src[i].coord[1] - src->low[1]) 
                         / (src->high[1] - src->low[1]));
          if (yi==ny) yi = ny-1;
          zi = (int) (nz * (params->src[i].coord[2] - src->low[2]) 
                         / (src->high[2] - src->low[2]));
          if (zi==nz) zi = nz-1;

          for (ti = 0; ti<nt; ti++)
          { Q = params->src[i].Q;
            if (nt>1)
            { double l, h, d;
              l = ti * highest_time / nt;
              h = (ti+1) * highest_time / nt;
              d = 0;
              if (params->src[i].start>l)
              { d += params->src[i].start-l;
              }
              if (params->src[i].stop<h)
              { d += h-params->src[i].stop;
              }
              if (d>h-l) continue;
              Q *= 1 - d/(h-l);
            }
            cells [xi + nx * (yi + ny * (zi + nz*ti))] += f*Q;
          }
        }

        N_records_used += 1;
      }

      /* Find next target index, if going by number of records. */

      if (mod_no<0) target_index += (double) (hindex-lindex) / (-mod_no-1);
    }

    log_file_close(&logfile);

    /* See if we're done. */

    if (*ap==0) break;

    /* Open next log file.  Note that initial records (architecture, etc.) 
       are read only from the first log file. */

    logfile.file_name = *ap++;
    log_file_open(&logfile,0);

    logg.index['q'] = -1;  /* So it won't be mistaken for a record from */
                           /* the new log file.                         */
  }

  if (N_records_used==0) 
  { fprintf(stderr,"None of the specified iterations were found\n");
    exit(1);
  }

  /* Output table of cell amounts. */

  for (xi = 0; xi<nx; xi++)
  { for (yi = 0; yi<ny; yi++)
    { for (zi = 0; zi<nz; zi++)
      { for (ti = 0; ti<nt; ti++)
        { printf("%12.4f", 
            cells [xi + nx * (yi + ny * (zi + nz*ti))] / sum_weights);
          if (hx) 
          { printf (" %+12.4f",
                    src->low[0]+(xi+0.5)*(src->high[0]-src->low[0])/nx);
          }
          if (hy) 
          { printf (" %+12.4f",
                    src->low[1]+(yi+0.5)*(src->high[1]-src->low[1])/ny);
          }
          if (hz) 
          { printf (" %+12.4f",
                    src->low[2]+(zi+0.5)*(src->high[2]-src->low[2])/nz);
          }
          if (ht) 
          { printf (" %12.4f", (ti+0.5)*src->max_stop/nt);
          }
          printf("\n");
        }
      }
    }
  }

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf (stderr, 
           "src-intensity n-x [ n-y [ n-z [ n-t ] ] ] / { log-file range }\n");
  exit(1);
}
