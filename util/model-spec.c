/* MODEL-SPEC.C - Program to specify a data model. */

/* Copyright (c) 1996 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "log.h"
#include "prior.h"
#include "model.h"


static void usage(void);


/* MAIN PROGRAM. */

void main
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

      case 'C': 
      { printf("class");   
        break;
      }

      case 'R': 
      { printf("real %s", prior_show(ps,m->noise));
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
              { printf("\n  %2d %9.4lf",i+1,v->time[i]);
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
   "Usage: model-spec log-file ( real noise-prior | binary | class | survival ... )\n");

  fprintf(stderr,
   "   or: model-spec log-file (to display stored specifications)\n");

  exit(1);
}

