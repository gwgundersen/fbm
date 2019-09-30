/* MODEL-SPEC.C - Program to specify a data model. */

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

#include "log.h"
#include "prior.h"
#include "model.h"
#include "matrix.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  static model_specification model, *m = &model; /* Static so that irrelevant */
  static model_survival      surv,  *v = &surv;  /*   fields are set to zero  */

  log_file logf;
  log_gobbled logg;

  char ps[100];
  char **ap;
  int i;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing specifications. */

  if (argc==2)
  {
    /* Open log file and gobble up initial records. */
  
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);

    logg.req_size['M'] = sizeof (model_specification);
    logg.req_size['V'] = sizeof (model_survival);

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }
  
    /* Display information in model record. */
  
    if ((m = logg.data['M'])==0)
    { printf ("\nNo model specification found\n\n");
      exit(0);
    }
  
    printf("\nData model:\n\n  ");
  
    switch (m->type)
    { 
      case 'B': 
      { printf("binary"); 
        break;
      }

      case 'N': 
      { printf("count"); 
        break;
      }

      case 'C': 
      { printf("class");   
        break;
      }

      case 'R': 
      { printf("real %s", prior_show(ps,m->noise));
        if (m->autocorr==1)
        { int i;
          printf(" acf");
          for (i = 0; i<m->n_autocorr; i++) printf(" %f",m->acf[i]);
        }
        break;
      }

      case 'V': 
      { 
        printf("survival "); 
 
        if ((v = logg.data['V'])==0)
        { printf("(no hazard type recorded)");
        }

        else
        { switch (v->hazard_type)
          { 
            case 'C': 
            { printf("const-hazard"); 
              break;
            }

            case 'P':
            { printf("pw-const-hazard");
              if (v->log_time) printf(" log");
              printf("\n\n");
              printf("  Time points for piecewise constant hazard:\n"); 
              for (i = 0; i<Max_time_points && v->time[i]>0; i++)
              { printf("\n  %2d %9.4f",i+1,v->time[i]);
              }
              break;
            }

            default: 
            { printf("(unknown hazard type)");
              break;
            }
          }
        }
        break;
      }
  
      default:  
      { printf("unknown type: %c",m->type); 
        break;
      }
    }

    printf("\n\n");
  
    log_file_close(&logf);
  
    exit(0);
  }

  /* Otherwise, examine arguments describing model. */

  m->type = 0;
  
  ap = argv+2;

  if (*ap==0) usage();
  
  if (strcmp(*ap,"binary")==0) 
  { m->type = 'B';
  }
  else if (strcmp(*ap,"count")==0) 
  { m->type = 'N';
  }
  else if (strcmp(*ap,"class")==0) 
  { m->type = 'C';
  }
  else if (strcmp(*ap,"real")==0)
  { m->type = 'R';
    if (*++ap==0 || !prior_parse(&m->noise,*ap)) usage();
    if (m->noise.scale || m->noise.two_point)
    { fprintf(stderr,"Illegal prior for noise level\n");
      exit(1);
    }
    ap += 1;
    if (*ap!=0)
    { if (strcmp(*ap,"acf")==0)
      { m->autocorr = 1;
        if (m->noise.alpha[2]!=0)
        { fprintf(stderr,
     "Autocorrelated noise is not allowed with case-by-case noise variances\n");
          exit(1);
        }
        ap += 1;
        if (*ap==0) usage();
        m->n_autocorr = 0;
        while (*ap!=0)
        { if (m->n_autocorr>=Max_autocorr) 
          { fprintf(stderr,
              "Autocorrelations are specified to too high a lag (max %d)\n",
              Max_autocorr);
            exit(1);
          }
          m->acf[m->n_autocorr] = atof(*ap);
          if (m->acf[m->n_autocorr]<-1 || m->acf[m->n_autocorr]>1
           || m->acf[m->n_autocorr]==0 && **ap!='0')
          { fprintf(stderr, "Invalid autocorrelation specified\n");
            exit(1);
          }
          m->n_autocorr += 1;
          ap += 1;
        }
        /* Check that the autocorrelation function is positive definite. This
           is done by trying it out on a series of length 100, though this
           is not completely guaranteed to find the problem if there is one. */
        { double cv[100][100];
          int i, j;
          for (i = 0; i<100; i++)
          { for (j = 0; j<i; j++)
            { cv[i][j] = i-j>m->n_autocorr ? 0 : m->acf[i-j-1];
            }
            cv[i][i] = 1;
          }
          if (!cholesky(&cv[0][0],100,0))
          { fprintf (stderr,
              "The specified autocorrelations are not positive definite\n");
            exit(1);
          }
        }
      }
      else
      { usage();
      }
    }
    ap -= 1;
  }
  else if (strcmp(*ap,"survival")==0) 
  { 
    m->type = 'V';
    ap += 1;

    if (*ap==0) usage();

    if (strcmp(*ap,"const-hazard")==0)
    { v->hazard_type = 'C';
    }

    else if (strcmp(*ap,"pw-const-hazard")==0)
    { 
      v->hazard_type = 'P';
      v->log_time = 0;
      if (*(ap+1)!=0 && strcmp(*(ap+1),"log")==0)
      { v->log_time = 1;
        ap += 1;
      }

      for (i = 0; *(ap+1)!=0; i++)
      { if (i==Max_time_points)
        { fprintf (stderr,
            "Too many time points with pw-const-hazard (max %d)\n",
            Max_time_points);
          exit(1);
        }
        v->time[i] = atof(*(ap+1));
        if (v->time[i]<=0 || i>0 && v->time[i]<=v->time[i-1])
        { fprintf(stderr,"Bad time point for pw-const-hazard: %s\n",*(ap+1));
          exit(1);
        }
        ap += 1;
      }

      v->time[i] = 0;
      if (i<2)
      { fprintf(stderr,
          "At least two time points needed with pw-const-hazard\n");
        exit(1);
      }
    }

    else 
    { fprintf(stderr,"Unknown hazard type: %s\n",*ap);
      exit(1);
    }
  }
  else
  { fprintf(stderr,"Unknown data model: %s\n",*ap);
    exit(1);
  }

  ap += 1;

  if (*ap!=0) usage();

  /* Append model (and maybe survival) record to log file. */

  log_file_open (&logf, 1);

  logf.header.type = 'M';
  logf.header.index = -1;
  logf.header.size = sizeof *m;
  log_file_append(&logf,m);

  if (m->type=='V')
  { logf.header.type = 'V';
    logf.header.index = -1;
    logf.header.size = sizeof *v;
    log_file_append(&logf,v);
  }

  log_file_close(&logf);

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
   "Usage: model-spec log-file model-specification...\n");
  fprintf(stderr,
   "   or: model-spec log-file (to display stored specifications)\n");
  fprintf(stderr, "Model specification:\n");
  fprintf(stderr, "   real noise-prior [ \"acf\" corr { corr } ]\n");
  fprintf(stderr, "      | binary | count | class | survival ... \n");

  exit(1);
}
