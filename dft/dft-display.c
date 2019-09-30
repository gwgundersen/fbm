/* DFT-DISPLAY.C - Program to print parameters, etc. of diffusion tree model. */

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
#include "dft.h"


static void usage(void)
{ fprintf (stderr, 
   "dft-display [ -b ] [ -h ] [ -p ] [ -l ] [ -t ] [ -n | -N ] [ -g[time] ]\n");
  fprintf (stderr,
   "            log-file [ index ]\n");
  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  dft_spec *dft;
  dft_hypers *h;
  model_specification *m;

  log_file logf;
  log_gobbled logg;

  int show_hypers;
  int show_params;
  int show_latent;
  int show_trees;
  int show_nodes;
  int show_graph;
  int include_locations;
  double graph_time;
  int bare;

  int index;

  int *parents;	
  double *divt;	
  double *latent;
  double *locations;

  int N_train;

  dft_state st;
  dft_tree_node *nodes;

  int dt;

  /* Look at arguments. */

  show_hypers= 0;
  show_params = 0;
  show_latent = 0;
  show_trees = 0;
  show_graph = 0;
  show_nodes = 0;
  include_locations = 0;
  bare = 0;

  while (argc>1 && *argv[1]=='-')
  {
    if (strcmp(argv[1],"-h")==0)
    { show_hypers = 1;
    }
    else if (strcmp(argv[1],"-p")==0)
    { show_params = 1;
    }
    else if (strcmp(argv[1],"-l")==0)
    { show_latent = 1;
    }
    else if (strcmp(argv[1],"-t")==0)
    { show_trees = 1;
    }
    else if (strcmp(argv[1],"-n")==0)
    { if (show_nodes) usage();
      show_nodes = 1;
      include_locations = 0;
    }
    else if (strcmp(argv[1],"-N")==0)
    { if (show_nodes) usage();
      show_nodes = 1;
      include_locations = 1;
    }
    else if (strcmp(argv[1],"-b")==0)
    { bare = 1;
    }
    else if (argv[1][0]=='-' && argv[1][1]=='g')
    { double gt;
      gt = argv[1][2]==0 ? 1.0 : atof(argv[1]+2);
      if (gt<=0 || gt>1) 
      { fprintf(stderr,"Time after -g must be in (0,1]\n");
        exit(1);
      }
      if (show_graph && gt!=graph_time)
      { fprintf(stderr,"More than one -g given with different times\n");
        exit(1);
      }
      graph_time = gt;
      show_graph = 1;
    }
    else 
    { usage();
    }

    argv += 1;
    argc -= 1;
  }

  if (!show_hypers && !show_params && !show_latent 
   && !show_trees && !show_graph && !show_nodes)
  { show_hypers = 1;
    show_params = 1;
  }

  index = -1;

  if (argc!=2 && argc!=3 
   || argc>2 && (index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0) 
  { usage();
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  dft_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  dft = logg.data['P'];
  m = logg.data['M'];

  dft_check_specs_present (dft, 0, m);

  /* Read the desired records from the log file. */

  if (index<0)
  { 
    log_gobble_last(&logf,&logg);

    if (logg.last_index<0)
    { fprintf(stderr,"No diffusion tree model record in log file\n");
      exit(1);
    }

    index = logg.last_index;
  }
  else
  {
    while (!logf.at_end && logf.header.index!=index)
    { log_file_forward(&logf);
    }

    if (logf.at_end)
    { fprintf(stderr,"That index does not appear in the log file\n");
      exit(1);
    }

    log_gobble(&logf,&logg);
  }

  if (logg.index['S']!=index)
  { fprintf(stderr,"No hyperparameters stored with that index\n");
    exit(1);
  }

  h = logg.data['S'];

  N_train = logg.data['T']==0 || logg.index['T']!=index ? 0 
          : logg.actual_size['T']/dft->N_trees/sizeof(double) + 1;

  dft_check_sizes (&logg, dft, m, N_train);

  /* Set up state. */

  parents   = logg.data['R']==0 || logg.index['R']!=index ? 0 : logg.data['R'];
  divt      = logg.data['T']==0 || logg.index['T']!=index ? 0 : logg.data['T'];
  latent    = logg.data['L']==0 || logg.index['L']!=index ? 0 : logg.data['L'];
  locations = logg.data['N']==0 || logg.index['N']!=index ? 0 : logg.data['N'];

  if (!parents)
  { nodes = 0;
  }
  else
  { nodes = chk_alloc (dft->N_trees*N_train, sizeof *nodes);

    for (dt = 0; dt<dft->N_trees; dt++)
    { (void) dft_conv_tree (parents+dt*(2*N_train)+N_train-1, 
                            nodes+dt*N_train, N_train);
    }
  }

  dft_setup_state (dft, m, st, h, parents, divt, locations, nodes, N_train);

  /* Print things that were asked for. */

  if (!bare) 
  { printf("\nDIFFUSION TREE MODEL IN FILE \"%s\" WITH INDEX %d\n\n", 
           logf.file_name, index);
  }

  if (show_hypers) 
  { dft_print_hypers (dft, m, h);
  }

  if (show_params)
  { dft_print_params (dft, h);
  }

  if (show_latent)
  { if (latent==0)
    { if (!bare) printf("\nNO LATENT VECTORS STORED WITH THIS INDEX\n\n");
    }
    else
    { dft_print_latent (dft, latent, N_train, bare);
    }
  }

  if (show_trees)
  { if (parents==0 || divt==0)
    { if (!bare) printf("\nNO DIFFUSION TREES STORED WITH THIS INDEX\n\n");
    }
    else
    { dft_print_trees (dft, st, N_train, bare);
    }
  }

  if (show_nodes)
  { if (parents==0 || divt==0)
    { if (!bare && !show_trees) 
      { printf("\nNO DIFFUSION TREES STORED WITH THIS INDEX\n\n");
      }
    }
    else
    { dft_print_nodes (dft, st, N_train, include_locations, bare);
    }
  }

  if (show_graph)
  { if (parents==0 || divt==0)
    { if (!bare) printf("\nNO DIFFUSION TREES STORED WITH THIS INDEX\n\n");
    }
    else
    { dft_graph_trees (dft, st, N_train, graph_time, bare);
    }
  }

  exit(0);
}
