/* NET-EVAL.C - Program to evaluate the network function at a grid of points. */

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
#include "data.h"
#include "model.h"
#include "net.h"

static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  net_arch *a;
  net_flags *flgs;
  model_specification *m;
  model_survival *sv;

  net_sigmas sigmas, *s = &sigmas;
  net_params params, *w = &params;
  net_values values, *v = &values;

  net_value *value_block;
  int value_count;

  static net_value grid_low[Max_inputs];
  static net_value grid_high[Max_inputs];

  static int grid_size[Max_inputs];
  static int grid_point[Max_inputs];

  log_file logf;
  log_gobbled logg;

  int lindex, hindex, index_mod;
  int ng;

  double targets[Max_targets];
  int gen_targets;
  int N_targets;

  char **ap;
  int first;
  int i, j;

  /* Look at arguments. */
  
  gen_targets = 0;

  if (argc>1 && strcmp(argv[argc-1],"targets")==0)
  { gen_targets = 1;
    argc -= 1;
    argv[argc] = 0;
  }

  if (argc<7 || (argc-3)%4!=0 || (argc-3)/4>Max_inputs)
  { usage();
  }

  logf.file_name = argv[1];
  
  parse_range (argv[2], &lindex, &hindex, &index_mod);

  if (hindex==-2) 
  { hindex = lindex;
  }
  else if (hindex==-1)
  { hindex = 1000000000;
  }

  ng = 0;

  for (ap = argv+3; *ap!=0; ap += 4)
  { if (strcmp(ap[0],"/")!=0 
     || ((grid_size[ng] = atoi(ap[3]))<=0 && strcmp(ap[3],"0")!=0)) usage();
    grid_low[ng] = atof(ap[1]);
    grid_high[ng] = atof(ap[2]);
    ng += 1;
  }

  /* Open log file and read network architecture and data model. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  net_record_sizes(&logg);

  if (!logf.at_end && logf.header.index==-1)
  { log_gobble(&logf,&logg);
  }

  a = logg.data['A'];
  flgs = logg.data['F'];
  m = logg.data['M'];
  sv = logg.data['V'];

  if (a==0)
  { fprintf(stderr,"No architecture specification in log file\n");
    exit(1);
  }

  if (gen_targets) 
  { 
    if (m==0)
    { fprintf(stderr,"No model specification in log file\n");
      exit(1);
    }

    if (m->type=='V' && v==0)
    { fprintf(stderr,"No hazard specification for survival model\n");
      exit(1);
    }

    if (m->type=='V' && sv->hazard_type!='C')
    { fprintf(stderr,
 "Can't generate targets randomly for survival model with non-constant hazard\n"
      );
      exit(1);
    }

    N_targets = model_targets(m,a->N_outputs);
  }

  if (a->N_inputs!=ng)
  { fprintf(stderr,
      "Number of grid ranges doesn't match number of input dimensions\n");
    exit(1);
  }

  s->total_sigmas = net_setup_sigma_count(a,flgs,m);
  w->total_params = net_setup_param_count(a,flgs);

  logg.req_size['S'] = s->total_sigmas * sizeof(net_sigma);
  logg.req_size['W'] = w->total_params * sizeof(net_param);

  /* Allocate space for values in network. */

  value_count = net_setup_value_count(a);
  value_block = chk_alloc (value_count, sizeof *value_block);

  net_setup_value_pointers (v, value_block, a);

  /* Evaluate function for the specified networks. */

  first = 1;

  for (;;)
  {
    /* Skip to next desired index, or to end of range. */
    
    while (!logf.at_end && logf.header.index<=hindex
     && (logf.header.index<lindex 
          || (logf.header.index-lindex)%index_mod!=0))
    { log_file_forward(&logf);
    }
  
    if (logf.at_end || logf.header.index>hindex)
    { break;
    }
  
    /* Gobble up records for this index. */

    log_gobble(&logf,&logg);
       
    /* Read the desired network from the log file. */
  
    if (logg.data['W']==0 || logg.index['W']!=logg.last_index)
    { fprintf(stderr,
        "No weights stored for the network with index %d\n",logg.last_index);
      exit(1);
    }
  
    w->param_block = logg.data['W'];
    net_setup_param_pointers (w, a, flgs);
  
    if (gen_targets) 
    {
      if (logg.data['S']==0 || logg.index['S']!=logg.last_index)
      { fprintf(stderr,
          "No sigmas stored for network with index %d\n",logg.last_index);
        exit(1);
      }
  
      s->sigma_block = logg.data['S'];
      net_setup_sigma_pointers (s, a, flgs, m);
    }
  
    /* Print the value of the network function, or targets generated from it, 
       for the grid of points. */

    if (first)
    { first = 0;
    }
    else
    { printf("\n");
    }
  
    for (i = 0; i<a->N_inputs; i++) 
    { grid_point[i] = 0;
      v->i[i] = grid_low[i];
    }
  
    for (;;)
    {
      net_func (v, 0, a, flgs, w);
  
      for (i = 0; i<a->N_inputs; i++) printf(" %8.5f",v->i[i]);
  
      if (gen_targets)
      { net_model_guess (v, targets, a, flgs, m, sv, w, s, 1);
        for (j = 0; j<N_targets; j++) printf(" %+.6e",targets[j]);
      }
      else
      { for (j = 0; j<a->N_outputs; j++) printf(" %+.6e",v->o[j]);
      }
  
      printf("\n");
  
      for (i = a->N_inputs-1; i>=0 && grid_point[i]==grid_size[i]; i--) 
      { grid_point[i] = 0;
        v->i[i] = grid_low[i];
      }
   
      if (i<0) break;
  
      if (i!=a->N_inputs-1) printf("\n");
  
      grid_point[i] += 1;
      v->i[i] = grid_low[i] 
                  + grid_point[i] * (grid_high[i]-grid_low[i]) / grid_size[i];
    }
  }
  
  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf (stderr, 
  "Usage: net-eval log-file range { / low high grid-size } [ \"targets\" ]\n");
  exit(1);
}
