/* NET-SPEC.C - Program to specify a new network (and create log file). */

/* Copyright (c) 1995, 1996 by Radford M. Neal 
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

#include "misc.h"
#include "log.h"
#include "prior.h"
#include "model.h"
#include "net.h"


static void usage(void);


/* MAIN PROGRAM. */

void main
( int argc,
  char **argv
)
{
  static net_arch    arch,   *a = &arch;   /* Static so irrelevant fields are */
  static net_priors  priors, *p = &priors; /*   set to zero (just in case)    */

  log_file logf;
  log_gobbled logg;

  char ps[100];
  char **ap;
  int l;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display specifications for existing network. */

  if (argc==2)
  {
    /* Open log file and gobble up initial records. */
  
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    net_record_sizes(&logg);

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }
  
    /* Display architecture record. */
  
    printf("\n");
  
    if ((a = logg.data['A'])==0)
    { printf ("No architecture specification found\n\n");
      exit(0);
    }

    printf("Network Architecture:\n\n");
  
    printf ("  Size of input layer:    %d\n", a->N_inputs);  
    if (a->N_layers>0)
    { printf ("  Sizes of hidden layers:");
      for (l = 0; l<a->N_layers; l++) 
      { printf(" %d",a->N_hidden[l]);
      }
      printf("\n");
    }
    printf ("  Size of output layer:   %d\n", a->N_outputs);
  
    printf("\n");
    
    /* Display priors record. */
  
    printf("\n");
  
    if ((p = logg.data['P'])==0)
    { printf("No prior specifications found\n\n");
      exit(0);
    }
  
    printf("Prior Specifications:\n");
  
    if (a->has_ti) 
    { printf("\n  Input Offsets:          %s\n",prior_show(ps,p->ti));
    }
  
    for (l = 0; l<a->N_layers; l++)
    { printf("\n         Hidden Layer %d\n\n",l);
      if (l>0 && a->has_hh[l-1])
      { printf("  Hidden-Hidden Weights:  %s\n", prior_show(ps,p->hh[l-1]));
      }
      if (a->has_ih[l]) 
      { printf("  Input-Hidden Weights:   %s\n", prior_show(ps,p->ih[l]));
      }
      if (a->has_bh[l]) 
      { printf("  Hidden Biases:          %s\n", prior_show(ps,p->bh[l]));
      }
      if (a->has_th[l]) 
      { printf("  Hidden Offsets:         %s\n", prior_show(ps,p->th[l]));
      }
    }

    printf("\n         Output Layer\n\n");

    for (l = a->N_layers-1; l>=0; l--)
    { if (a->has_ho[l]) 
      { printf("  Hidden%d-Output Weights: %s\n",l,prior_show(ps,p->ho[l]));
      }
    }

    if (a->has_io) 
    { printf("  Input-Output Weights:   %s\n",prior_show(ps,p->io));
    }

    if (a->has_bo) 
    { printf("  Output Biases:          %s\n",prior_show(ps,p->bo));
    }

    for (l = 0; l<a->N_layers && !a->has_ah[l]; l++) ;

    if (l<a->N_layers || a->has_ao)
    {
       if (a->N_layers>0 && l<a->N_layers)
      { printf("\n  Hidden adjustments: ");
        for (l = 0; l<a->N_layers; l++)
        { if (p->ah[l]==0) printf(" -");
          else             printf(" %.2lf",p->ah[l]);
        }
      }

      if (a->has_ao)
      { printf("\n  Output adjustments: ");
        if (p->ao==0) printf(" -");
        else          printf(" %.2lf",p->ao);
      }
      printf("\n");
    }
  
    printf("\n");
  
    log_file_close(&logf);
  
    exit(0);
  }

  /* Otherwise, figure out architecture and priors from program arguments. */

  a->N_layers = 0;
  
  ap = argv+2;

  if (*ap==0 || (a->N_inputs = atoi(*ap++))<=0) usage();

  while (*ap!=0 && *(ap+1)!=0 && strcmp(*(ap+1),"/")!=0)
  { if (a->N_layers==Max_layers)
    { fprintf(stderr,"Too many layers specified (maximum is %d)\n",Max_layers);
      exit(1);
    }
    if ((a->N_hidden[a->N_layers] = atoi(*ap++))<=0) usage();
    a->N_layers += 1;
  }

  if (*ap==0 || (a->N_outputs = atoi(*ap++))<=0) usage();

  if (*ap==0 || strcmp(*ap,"/")!=0) usage();

  if (*++ap==0 || (a->has_ti = strcmp(*ap,"-")!=0) 
                    && !prior_parse(&p->ti,*ap)) usage();

  if (a->N_layers>0)
  { 
    for (l = 0; l<a->N_layers; l++)
    { 
      if (l>0)
      { if (*++ap==0 || (a->has_hh[l-1] = strcmp(*ap,"-")!=0)
                          && !prior_parse(&p->hh[l-1],*ap)) usage();
      }

      if (*++ap==0 || (a->has_ih[l] = strcmp(*ap,"-")!=0)
                        && !prior_parse(&p->ih[l],*ap)) usage();

      if (*++ap==0 || (a->has_bh[l] = strcmp(*ap,"-")!=0)
                        && !prior_parse(&p->bh[l],*ap)) usage();

      if (*++ap==0 || (a->has_th[l] = strcmp(*ap,"-")!=0)
                        && !prior_parse(&p->th[l],*ap)) usage();

    }

    for (l = a->N_layers-1; l>=0; l--)
    { if (*(ap+1)==0 || strcmp(*(ap+1),"/")==0
       || *(ap+2)==0 || strcmp(*(ap+2),"/")==0
       || *(ap+3)==0 || strcmp(*(ap+3),"/")==0)
      { a->has_ho[l] = 0;
      }
      else
      { if (*++ap==0 || (a->has_ho[l] = strcmp(*ap,"-")!=0)
                          && !prior_parse(&p->ho[l],*ap)) usage();
      }
    }
  }

  if (*++ap==0 || (a->has_io = strcmp(*ap,"-")!=0)
                    && !prior_parse(&p->io,*ap)) usage();

  if (*++ap==0 || (a->has_bo = strcmp(*ap,"-")!=0)
                    && !prior_parse(&p->bo,*ap)) usage();

  if (*++ap!=0 && strcmp(*ap,"/")==0 
   && *(ap+1)!=0 && strchr("+-.0123456789",**(ap+1))!=0)
  { 
    ap += 1;

    for (l = 0; l<a->N_layers; l++)
    { if (*ap==0) usage();
      p->ah[l] = 0;
      a->has_ah[l] = strcmp(*ap,"-")!=0;
      if (a->has_ah[l] && (p->ah[l] = atof(*ap))<=0) usage();
      ap += 1;
    }

    if (*ap==0) usage();
    p->ao = 0;
    a->has_ao = strcmp(*ap,"-")!=0;
    if (a->has_ao && (p->ao = atof(*ap))<=0) usage();
    ap += 1;
  }

  if (*ap!=0) usage();

  if (p->ti.scale || p->ti.alpha[2]!=0
   || p->bo.scale || p->bo.alpha[2]!=0)
  { fprintf(stderr,"Illegal prior for biases or offsets\n");
    exit(1); 
  }

  for (l = 0; l<a->N_layers; l++)
  { if (p->bh[l].scale || p->bh[l].alpha[2]!=0
     || p->th[l].scale || p->th[l].alpha[2]!=0)
    { fprintf(stderr,"Illegal prior for biases or offsets\n");
      exit(1); 
    }
  }

  /* Create log file and write records. */

  log_file_create(&logf);

  logf.header.type = 'A';
  logf.header.index = -1;
  logf.header.size = sizeof *a;
  log_file_append(&logf,a);

  logf.header.type = 'P';
  logf.header.index = -1;
  logf.header.size = sizeof *p;
  log_file_append(&logf,p);

  log_file_close(&logf);

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
   "Usage: net-spec log-file N-inputs { N-hidden } N-outputs \n");

  fprintf(stderr,
   "                / ti [ ih bh th { hh ih bh th } ] { ho } io bo  [ / { ah } ao ]\n");

  fprintf(stderr,
   "   or: net-spec log-file (to display stored network specifications)\n");

  fprintf(stderr,
   "Prior: [x]Width[:[Alpha-type][:[Alpha-unit][:[Alpha-weight]]]]\n");

  exit(1);
}

