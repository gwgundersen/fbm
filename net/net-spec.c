/* NET-SPEC.C - Program to specify a new network (and create log file). */

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
#include "prior.h"
#include "model.h"
#include "data.h"
#include "net.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  static net_arch    arch,   *a = &arch;   /* Static so irrelevant fields are */
  static net_priors  priors, *p = &priors; /*   set to zero (just in case)    */
  static net_flags   flags,  *flgs = &flags;
  int any_flags;

  log_file logf;
  log_gobbled logg;

  char ps[1000];
  char **ap;
  int i, j, l;

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

    flgs = logg.data['F'];

    printf("Network Architecture:\n\n");
  
    printf ("  Input layer:     size %d\n", a->N_inputs);  

    for (l = 0; l<a->N_layers; l++) 
    { printf("  Hidden layer %d:  size %d",l,a->N_hidden[l]);
      if (flgs==0 || flgs->layer_type[l]==Tanh_type) printf("  tanh");
      else if (flgs->layer_type[l]==Identity_type)   printf("  identity");
      else if (flgs->layer_type[l]==Sin_type)        printf("  sin");
      else                                           printf("  UNKNOWN TYPE!");
      if (flgs!=0 && l<7 
       && list_flags (flgs->omit, a->N_inputs, 1<<(l+1), ps) > 0)
      { printf("  omit%s",ps);
      }
      printf("\n");
    }

    printf ("  Output layer:    size %d", a->N_outputs);
    if (flgs!=0 && l<7 
     && list_flags (flgs->omit, a->N_inputs, 1, ps) > 0)
    { printf("  omit%s",ps);
    }
    printf("\n");
  
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
          else             printf(" %.2f",p->ah[l]);
        }
      }

      if (a->has_ao)
      { printf("\n  Output adjustments: ");
        if (p->ao==0) printf(" -");
        else          printf(" %.2f",p->ao);
      }
      printf("\n");
    }
  
    printf("\n");
  
    log_file_close(&logf);
  
    exit(0);
  }

  /* Otherwise, figure out architecture and priors from program arguments. */

  any_flags = 0;

  a->N_layers = 0;
  
  ap = argv+2;

  if (*ap==0 || (a->N_inputs = atoi(*ap++))<=0) usage();

  while (*ap!=0 && strcmp(*ap,"/")!=0)
  { 
    double size;
    int omit, type;
    int i;

    if ((size = atoi(*ap++))<=0) usage();

    if (*ap==0) usage();

    omit = 0;
    type = -1;

    while ((*ap)[0]>='a' && (*ap)[0]<='z')
    { if ((*ap)[0]=='o' && (*ap)[1]=='m' && (*ap)[2]=='i' && (*ap)[3]=='t' 
       && (*ap)[4]==':')
      { if (omit) usage();
        omit = 1;
        parse_flags (*ap+4, flgs->omit, a->N_inputs, 1);
      }
      else if (strcmp(*ap,"tanh")==0)
      { if (type>=0) usage();
        type = Tanh_type;
      }
      else if (strcmp(*ap,"identity")==0)
      { if (type>=0) usage();
        type = Identity_type;
      }
      else if (strcmp(*ap,"sin")==0)
      { if (type>=0) usage();
        type = Sin_type;
      }
      else
      { usage();
      }
      any_flags = 1;
      ap += 1;
    }

    if (*ap!=0 && strcmp(*ap,"/")!=0)
    { if (a->N_layers == (Max_layers>7 ? 7 : Max_layers))
      { fprintf(stderr,"Too many layers specified (maximum is %d)\n",
                        Max_layers>7 ? 7 : Max_layers);
        exit(1);
      }
      a->N_hidden[a->N_layers] = size;
      flgs->layer_type[a->N_layers] = type==-1 ? Tanh_type : type;
      for (i = 0; i<a->N_inputs; i++) 
      { flgs->omit[i] = (flgs->omit[i] | ((flgs->omit[i]&1)<<(a->N_layers+1))) & ~1;
      }
      a->N_layers += 1;
    }
    else
    { a->N_outputs = size;
      if (type!=-1) usage();
    }
  }

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

  if (any_flags)
  { logf.header.type = 'F';
    logf.header.index = -1;
    logf.header.size = sizeof *flgs;
    log_file_append(&logf,flgs);
  }

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
   "Usage: net-spec log-file N-inputs { N-hidden { flag } } N-outputs { flag }\n");

  fprintf(stderr,
   "                / ti [ ih bh th { hh ih bh th } ] { ho } io bo  [ / { ah } ao ]\n");

  fprintf(stderr,
   "   or: net-spec log-file (to display stored network specifications)\n");

  fprintf(stderr,
   "Prior: [x]Width[:[Alpha-type][:[Alpha-unit][:[Alpha-weight]]]]\n");

  exit(1);
}

