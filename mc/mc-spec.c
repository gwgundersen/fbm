/* MC-SPEC.C - Specify parameters of Markov chain Monte Carlo simulation. */

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
#include "mc.h"


static void usage (void);
static void display_specs (log_gobbled *);


/* MAIN PROGRAM. */

void main
( int argc,
  char **argv
)
{
  static mc_ops ops0, *ops = &ops0; /* Static so types get set to zero (empty)*/
  static mc_traj traj0, *traj = &traj0;

  int have_traj, any_specs;
  int o;

  log_file logf;
  log_gobbled logg;
  int index;
  int depth;

  char **ap;
  char *s;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing specifications. */

  if (argc==2 || argc==3 && *argv[2]>='0' && *argv[2]<='9')
  {
    index = argc==2 ? -1 : atoi(argv[2]);

    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    logg.req_size['o'] = sizeof *ops;
    logg.req_size['t'] = sizeof *traj;

    any_specs = 0;

    if (index==-1)
    {
      while (!logf.at_end)
      { log_gobble(&logf,&logg);
      }
  
      if (logg.data['o']!=0 || logg.data['t']!=0)
      { printf("\nLAST SPECIFICATIONS IN LOG FILE:\n");
        display_specs(&logg);
        any_specs = 1;
      }
  
    }
    else
    {
      while (!logf.at_end && logf.header.index<=index)
      { log_gobble(&logf,&logg);
      }

      if (logg.index['o']!=index) logg.data['o'] = 0;
      if (logg.index['t']!=index) logg.data['t'] = 0;
  
      if (logg.data['o']!=0 || logg.data['t']!=0)
      { printf("\nSPECIFICATIONS AT INDEX %d IN LOG FILE:\n",index);
        display_specs(&logg);
        any_specs = 1;
      }
    }

    if (!any_specs)
    { printf("\nNo Monte Carlo specification found\n\n");
    }

    exit(0);
  }

  /* Otherwise, look at remaining arguments. */

  ap = argv+2; 

  o = 0;

  while (*ap && strcmp(*ap,"/")!=0)
  {
    if (o==Max_mc_ops)
    { fprintf(stderr,"Too many operations (max %d)\n",Max_mc_ops);
      exit(1);
    }

    if (strcmp(*ap,"end")==0)
    { 
      ops->op[o].type = 'E';
      
      ap += 1;

      depth -= 1;
      if (depth<0) 
      { fprintf(stderr,"Too many 'end' operations\n");
        exit(1);
      }
    }

    else if (strcmp(*ap,"heatbath")==0)
    { 
      ops->op[o].type = 'B';
      ops->op[o].heatbath_decay = 0;

      ap += 1;

      if (*ap && strchr("0123456789+-.",**ap))
      { if ((ops->op[o].heatbath_decay = atof(*ap++))<0) usage();
      }
    }

    else if (strcmp(*ap,"radial-heatbath")==0)
    { 
      ops->op[o].type = 'r';
      ap += 1;
    }

    else if (strcmp(*ap,"negate")==0)
    { 
      ops->op[o].type = 'N';
      ap += 1;
    }

    else if (strcmp(*ap,"metropolis")==0)
    {
      ops->op[o].type = 'M';
      ops->op[o].stepsize_adjust = 1;
      ops->op[o].stepsize_alpha = 0;

      ap += 1;

      if (*ap && strchr("0123456789+-.",**ap))
      { if ((ops->op[o].stepsize_adjust = atof(*ap))==0) usage();
        if (strchr(*ap,':')!=0)
        { if ((ops->op[o].stepsize_alpha = atof(strchr(*ap,':')+1))==0) usage();
        }
        ap += 1;
      }
    }

    else if (strcmp(*ap,"met-1")==0 
          || strcmp(*ap,"slice-1")==0
          || strcmp(*ap,"slice-over")==0)
    {
      ops->op[o].type = strcmp(*ap,"met-1")==0 ? 'm'
                      : strcmp(*ap,"slice-1")==0 ? 'S' : 'O';
                      
      ops->op[o].refinements = 0;
      ops->op[o].refresh_prob = 0.0;
      ops->op[o].stepsize_adjust = 1;
      ops->op[o].stepsize_alpha = 0;
      ops->op[o].steps = 0;
      ops->op[o].firsti = -1;

      ap += 1;

      if (ops->op[o].type=='O')
      { 
        if (*ap && strchr("0123456789",**ap))
        { if ((ops->op[o].refinements = atoi(*ap)) < 0) usage();
          ap += 1;
        }

        if (*ap && strchr("0123456789.",**ap))
        { if ((ops->op[o].refresh_prob = atof(*ap)) < 0) usage();
          ap += 1;
        }
      }

      if (*ap && strchr("0123456789+-.",**ap))
      { if ((ops->op[o].stepsize_adjust = atof(*ap))==0) usage();
        if (strchr(*ap,':')!=0)
        { if ((ops->op[o].stepsize_alpha = atof(strchr(*ap,':')+1))==0) usage();
        }
        ap += 1;
      }

      if (*ap && strchr("0123456789",**ap))
      { ops->op[o].steps = atoi(*ap);
        ap += 1;
      }

      if (*ap && strchr("0123456789",**ap))
      { ops->op[o].firsti = atoi(*ap);
        ops->op[o].lasti = ops->op[o].firsti;
        if (strchr(*ap,':')!=0)
        { if ((ops->op[o].lasti = atoi(strchr(*ap,':')+1)) < ops->op[o].firsti)
          { usage();
          }
        }
        ap += 1;
      }
    }

    else if (strcmp(*ap,"dynamic")==0 || strcmp(*ap,"permuted-dynamic")==0
          || strcmp(*ap,"slice-inside")==0 || strcmp(*ap,"slice-outside")==0)
    {
      ops->op[o].type = strcmp(*ap,"dynamic")==0 ? 'D' 
                      : strcmp(*ap,"permuted-dynamic")==0 ? 'P'
                      : strcmp(*ap,"slice-inside")==0 ? 'i' : 'o';

      ops->op[o].stepsize_adjust = 1;
      ops->op[o].stepsize_alpha = 0;

      ap += 1;

      if (!*ap || !strchr("0123456789+-.",**ap)) usage();

      if ((ops->op[o].steps = atoi(*ap))<=0) usage();

      if (ops->op[o].type=='o')
      { ops->op[o].in_steps = ops->op[o].steps;
        if (strchr(*ap,'/')!=0)
        { if ((ops->op[o].in_steps = atoi(strchr(*ap,'/')+1))<=0
           || ops->op[o].in_steps>ops->op[o].steps) usage();
        }
      }

      ap += 1;

      if (*ap && strchr("0123456789+-.",**ap))
      { if ((ops->op[o].stepsize_adjust = atof(*ap))==0) usage();
        if (strchr(*ap,':')!=0)
        { if ((ops->op[o].stepsize_alpha = atof(strchr(*ap,':')+1))==0) usage();
        }
        ap += 1;
      }
    }

    else if (strcmp(*ap,"hybrid")==0 || strcmp(*ap,"tempered-hybrid")==0)
    {
      if (strcmp(*ap,"hybrid")==0)
      {
        ops->op[o].type = 'H';
      }
      else
      {
        ops->op[o].type = 'T';

        ap += 1;
        if (!*ap || !strchr("0123456789+-.",**ap)) usage();
        
        if ((ops->op[o].temper_factor = atof(*ap))<=0) usage();
      }

      ops->op[o].stepsize_adjust = 1;
      ops->op[o].stepsize_alpha = 0;
      ops->op[o].in_steps = 0;
      ops->op[o].window = 1;
      ops->op[o].jump = 1;

      ap += 1;

      if (!*ap || !strchr("0123456789+-.",**ap)) usage();

      if ((ops->op[o].steps = atoi(*ap))<=0) usage();
      if (strchr(*ap,'/')!=0)
      { s = strchr(*ap,'/')+1;
        if ((ops->op[o].in_steps = atoi(s))<=0) usage();
        if (strchr(s,':')!=0)
        { if ((ops->op[o].jump = atoi(strchr(s,':')+1))<=0) usage();
          if (ops->op[o].steps%ops->op[o].jump!=0)
          { fprintf (stderr,
             "Maximum number of steps must be multiple of jump\n");
            exit(1);
          }
        }
      } 
      else if (strchr(*ap,':')!=0)
      { s = strchr(*ap,':')+1;
        if ((ops->op[o].window = atoi(s))<=0) usage();
        if (strchr(s,':')!=0)
        { if ((ops->op[o].jump = atoi(strchr(s,':')+1))<=0) usage();
          if (ops->op[o].steps%ops->op[o].jump!=0)
          { fprintf(stderr,"Total number of steps must be multiple of jump\n");
            exit(1);
          }
        }
        if (ops->op[o].window > ops->op[o].steps/ops->op[o].jump + 1)
        { fprintf(stderr,"Window can't be bigger than whole trajectory\n");
          exit(1);
        }
      } 
      ap += 1;

      if (*ap && strchr("0123456789+-.",**ap))
      { if ((ops->op[o].stepsize_adjust = atof(*ap))==0) usage();
        if (strchr(*ap,':')!=0)
        { if ((ops->op[o].stepsize_alpha = atof(strchr(*ap,':')+1))==0) usage();
        }
        ap += 1;
      }
    }

    else if (strcmp(*ap,"repeat")==0)
    {
      ops->op[o].type = 'R';

      ap += 1;

      if (!*ap || !strchr("0123456789+-.",**ap)) usage();

      if ((ops->op[o].repeat_count = atoi(*ap++))<=0) usage();

      depth += 1;
    }

    else if (strcmp(*ap,"temp-trans")==0)
    {
      ops->op[o].type = 't';

      ap += 1;

      depth += 1;
    }

    else if (strcmp(*ap,"sim-temp")==0)
    { ops->op[o].type = 's';
      ap += 1;
    }

    else if (strcmp(*ap,"rand-dir")==0)
    { 
      ops->op[o].type = 'b';
      ap += 1;
    }

    else if (strcmp(*ap,"neg-dir")==0)
    { 
      ops->op[o].type = 'n';
      ap += 1;
    }

    else if (strcmp(*ap,"plot")==0)
    { 
      ops->op[o].type = 'p';
      ap += 1;
    }

    else if (**ap>='a' && **ap<='z' || **ap>='A' && **ap<='Z')
    {
      ops->op[o].type = 'A';
      ops->op[o].app_param = 0;

      strcpy(ops->op[o].appl,*ap);

      ap += 1;

      if (*ap && strchr("0123456789+-.",**ap))
      { ops->op[o].app_param = atof(*ap);
        ap += 1;
      }
    }

    else
    { usage();
    }

    o += 1;
  }

  have_traj = 0;

  if (*ap)
  {
    have_traj = 1;
    ap += 1;
  
    if (*ap && strcmp(*ap,"leapfrog")==0)
    {
      traj->type = 'L';
      traj->halfp = 1;
      traj->N_approx = 1;

      ap += 1;

      if (strcmp(*ap,"halfp")==0)
      { traj->halfp = 1;
        ap += 1;
      }
      else if (strcmp(*ap,"halfq")==0)
      { traj->halfp = 0;
        ap += 1;
      }

      if (*ap && (traj->N_approx = atoi(*ap++))==0) usage();
    }

    else 
    { usage();
    }

    if (*ap) usage();
  }

  /* Open log file and read last record. */

  log_file_open (&logf, 1);
  log_file_last (&logf);

  index = logf.at_end || logf.header.index<0 ? 0 : logf.header.index;

  /* Write specifications to log file. */

  if (o>0)
  { 
    logf.header.type = 'o';
    logf.header.index = index;
    logf.header.size = sizeof *ops;
    log_file_append (&logf, ops);
  }

  if (have_traj)
  { 
    logf.header.type = 't';
    logf.header.index = index;
    logf.header.size = sizeof *traj;
    log_file_append (&logf, traj);
  }

  log_file_close (&logf);

  exit(0);
}


/* DISPLAY SPECIFICATIONS. */

static void display_specs
( log_gobbled *logg
)
{
  mc_ops *ops;
  mc_traj *traj;
  int depth;
  int o;

  if (logg->data['o']!=0)
  {
    printf("\nOperations making up each iteration:\n\n");

    ops = logg->data['o'];
    depth = 0;
  
    for (o = 0; o<Max_mc_ops && ops->op[o].type!=0; o++)
    {
      printf(" %*s", 2 * (depth - (ops->op[o].type=='E')), "");

      switch (ops->op[o].type)
      {
        case 'E':
        { printf(" end\n");
          depth -= 1;
          break;
        }

        case 'B':
        { printf(" heatbath");
          if (ops->op[o].heatbath_decay!=0)
          { printf(" %.4lf",ops->op[o].heatbath_decay);
          }
          printf("\n");
          break;
        }

        case 'r':
        { printf(" radial-heatbath\n");
          break;
        }

        case 'N':
        { printf(" negate\n");
          break;
        }
  
        case 'M':
        { printf(" metropolis");
          if (ops->op[o].stepsize_alpha!=0)
          { printf(" %.4lf:%.4lf",ops->op[o].stepsize_adjust,
                                  ops->op[o].stepsize_alpha);
          }
          else if (ops->op[o].stepsize_adjust!=1)
          { printf(" %.4lf",ops->op[o].stepsize_adjust);
          }
          printf("\n");
          break;
        }

        case 'm':
        { printf(" met-1");
          if (ops->op[o].stepsize_alpha!=0)
          { printf(" %.4lf:%.4lf",ops->op[o].stepsize_adjust,
                                  ops->op[o].stepsize_alpha);
          }
          else if (ops->op[o].stepsize_adjust!=1 || ops->op[o].firsti!=-1)
          { printf(" %.4lf",ops->op[o].stepsize_adjust);
          }
          if (ops->op[o].firsti!=-1)
          { printf(" %d",ops->op[o].firsti);
            if (ops->op[o].lasti!=ops->op[o].firsti)
            { printf(":%d",ops->op[o].lasti);
            }
          }
          printf("\n");
          break;
        }
  
        case 'D': case 'P': case 'i': case 'o':
        { printf (" %s", ops->op[o].type=='D' ? "dynamic" 
                       : ops->op[o].type=='P' ? "permuted-dynamic"
                       : ops->op[o].type=='i' ? "slice-inside" 
                                              : "slice-outside");
          printf (" %d", ops->op[o].steps);
          if (ops->op[o].type=='o' && ops->op[o].in_steps!=ops->op[o].steps)
          { printf ("/%d", ops->op[o].in_steps);
          }

          if (ops->op[o].stepsize_alpha!=0)
          { printf(" %.4lf:%.4lf",ops->op[o].stepsize_adjust,
                                  ops->op[o].stepsize_alpha);
          }
          else if (ops->op[o].stepsize_adjust!=1)
          { printf(" %.4lf",ops->op[o].stepsize_adjust);
          }
          printf("\n");
          break;
        }
  
        case 'H': 
        { printf(" hybrid %d",ops->op[o].steps);
          if (ops->op[o].in_steps!=0)
          { printf("/%d",ops->op[o].in_steps);
          }
          else
          { if (ops->op[o].window!=1 || ops->op[o].jump!=1)
            { printf(":%d",ops->op[o].window);
            }
          }
          if (ops->op[o].jump!=1)
          { printf(":%d",ops->op[o].jump);
          }
          if (ops->op[o].stepsize_alpha!=0)
          { printf(" %.4lf:%.4lf",ops->op[o].stepsize_adjust,
                                  ops->op[o].stepsize_alpha);
          }
          else if (ops->op[o].stepsize_adjust!=1)
          { printf(" %.4lf",ops->op[o].stepsize_adjust);
          }
          printf("\n");
          break;
        }

        case 'S':
        { printf(" slice-1");
          if (ops->op[o].stepsize_alpha!=0)
          { printf(" %.4lf:%.4lf",ops->op[o].stepsize_adjust,
                                  ops->op[o].stepsize_alpha);
          }
          else if (ops->op[o].stepsize_adjust!=1 
                || ops->op[o].steps!=0
                || ops->op[o].firsti!=-1)
          { printf(" %.4lf",ops->op[o].stepsize_adjust);
          }
          if (ops->op[o].steps!=0 || ops->op[o].firsti!=-1)
          { printf(" %d",ops->op[o].steps);
          }
          if (ops->op[o].firsti!=-1)
          { printf(" %d",ops->op[o].firsti);
            if (ops->op[o].lasti!=ops->op[o].firsti)
            { printf(":%d",ops->op[o].lasti);
            }
          }
          printf("\n");
          break;
        }
    
        case 'O':
        { printf(" slice-over");
          if (ops->op[o].refinements!=0
           || ops->op[o].refresh_prob!=0 
           || ops->op[o].stepsize_alpha!=0
           || ops->op[o].stepsize_adjust!=1 
           || ops->op[o].firsti!=-1)
          { printf(" %d",ops->op[o].refinements);
          }
          if (ops->op[o].refresh_prob!=0 
           || ops->op[o].stepsize_alpha!=0
           || ops->op[o].stepsize_adjust!=1 
           || ops->op[o].firsti!=-1)
          { printf(" %.5lf",ops->op[o].refresh_prob);
          }
          if (ops->op[o].stepsize_alpha!=0)
          { printf(" %.4lf:%.4lf",ops->op[o].stepsize_adjust,
                                  ops->op[o].stepsize_alpha);
          }
          else if (ops->op[o].stepsize_adjust!=1 
                || ops->op[o].steps!=0
                || ops->op[o].firsti!=-1)
          { printf(" %.4lf",ops->op[o].stepsize_adjust);
          }
          if (ops->op[o].steps!=0 || ops->op[o].firsti!=-1)
          { printf(" %d",ops->op[o].steps);
          }
          if (ops->op[o].firsti!=-1)
          { printf(" %d",ops->op[o].firsti);
            if (ops->op[o].lasti!=ops->op[o].firsti)
            { printf(":%d",ops->op[o].lasti);
            }
          }
          printf("\n");
          break;
        }
  
        case 'T': 
        { printf(" tempered-hybrid %.6lf %d",ops->op[o].temper_factor,
                                              ops->op[o].steps);
          if (ops->op[o].window!=1)
          { printf(":%d",ops->op[o].window);
            if (ops->op[o].jump!=0)
            { printf(":%d",ops->op[o].jump);
            }
          }
          if (ops->op[o].stepsize_alpha!=0)
          { printf(" %.4lf:%.4lf",ops->op[o].stepsize_adjust,
                                  ops->op[o].stepsize_alpha);
          }
          else if (ops->op[o].stepsize_adjust!=1)
          { printf(" %.4lf",ops->op[o].stepsize_adjust);
          }
          printf("\n");
          break;
        }

        case 's':
        { printf(" sim-temp\n");
          break;
        }

        case 'b':
        { printf(" rand-dir\n");
          break;
        }

        case 'n':
        { printf(" neg-dir\n");
          break;
        }
  
        case 'A':
        { printf(" %s",ops->op[o].appl);
          if (ops->op[o].app_param!=0) printf(" %.4lf",ops->op[o].app_param);
          printf("\n");
          break;
        }
  
        case 'R':
        { printf(" repeat %d",ops->op[o].repeat_count);
          printf("\n");
          depth += 1;
          break;
        }
  
        case 't':
        { printf(" temp-trans\n");
          depth += 1;
          break;
        }

        case 'p':
        { printf(" plot\n");
          break;
        }
  
        default:
        { printf(" unknown operation: %c\n",ops->op[o].type);
          break;
        }
      }
    }
  }

  if (logg->data['t']!=0)
  {
    printf("\nMethod for computing trajectories:\n\n");

    traj = logg->data['t'];
    
    switch (traj->type)
    {
      case 'L': 
      { printf ("  leapfrog %s %d\n", traj->halfp ? "halfp" : "halfq", 
                                      traj->N_approx);
        break;
      }

      default:
      { printf("  unknown method: %c\n",traj->type);
        break;
      }
    }
  }
 
  printf("\n");
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr, 
    "Usage: mc-spec log-file { operation-spec } [ / trajectory-spec ]\n");

  fprintf(stderr,
    "   or: mc-spec log-file [ index ] (to display stored specifications)\n");

  exit(1);
}
