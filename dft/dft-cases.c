/* DFT-CASES.C - Generate cases from a diffusion tree model. */

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
#include "rand.h"
#include "dft.h"
#include "dft-data.h"


static void usage(void)
{ fprintf (stderr,
"Usage: dft-cases [ -# ]\n");
  fprintf (stderr,
"         [ -b ] [ -h ] [ -p ] [ -l ] [ -t ] [ -n | -N ] [ -g[time] ]\n");
  fprintf (stderr,
"         log-file index output-file n-cases [ random-seed ] [ new-log-file ]\n"
);
  exit(1);
}


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  dft_spec *dft;
  dft_hypers *hyp;
  model_specification *m;

  int N_targets, N_trees;

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
  int only_tree;

  int have_trees;
  int have_latent;
  int have_locations;

  int *parents, *parents0;	
  double *divt, *divt0;	
  double *latent, *latent0;
  double *locations, *locations0;
  double *terminals[Max_trees];

  int index, n_cases, n_tot;

  double *noise_SD;

  dft_state st;
  dft_tree_node *nodes;

  char *new_log_file;
  char *output_file;
  FILE *output;

  int random_seed;

  int first[Max_trees];

  int N, i, j, j0, t, dt;

  /* Look at arguments. */

  new_log_file = 0;
  if (argc>0 && strchr("0123456789+-",argv[argc-1][0])==0)
  { new_log_file = argv[argc-1];
    argc -= 1;
  }

  show_hypers= 0;
  show_params = 0;
  show_latent = 0;
  show_trees = 0;
  show_graph = 0;
  show_nodes = 0;
  include_locations = 0;
  bare = 0;
  only_tree = 0;

  while (argc>1 && *argv[1]=='-')
  {
    if (argv[1][0]=='-' && atoi(argv[1]+1)>0)
    { if (only_tree!=0) usage();
      only_tree = atoi(argv[1]+1);
    }
    else if (strcmp(argv[1],"-h")==0)
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

  if (argc!=5 && argc!=6) usage();

  logf.file_name = argv[1];

  if ((index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0) 
  { usage();
  }

  output_file = argv[3];

  if ((n_cases = atoi(argv[4]))<=0) usage();

  random_seed = index;
  if (argc==6 && (random_seed = atoi(argv[5]))<=0) usage();

  /* Open log file and read specification. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  dft_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  dft = logg.data['P'];
  m = logg.data['M'];

  dft_check_specs_present(dft,0,m);

  N_targets = dft->N_targets;
  N_trees = dft->N_trees;

  if (N_trees<only_tree)
  { fprintf(stderr,"Number of tree selected for display is too big\n");
    exit(1);
  }

  /* Read the desired records from the log file. */

  while (!logf.at_end && logf.header.index!=index)
  { log_file_forward(&logf);
  }

  if (logf.at_end)
  { fprintf(stderr,"That index does not appear in the log file\n");
    exit(1);
  }

  log_gobble(&logf,&logg);

  if (logg.index['S']!=index)
  { fprintf(stderr,"No hyperparametes stored with that index\n");
    exit(1);
  }

  hyp = logg.data['S'];

  noise_SD = hyp->SD + N_trees*N_targets;

  have_trees     = logg.data['T']!=0;
  have_latent    = logg.data['L']!=0;
  have_locations = logg.data['N']!=0;

  if (have_trees && logg.data['R']==0)
  { fprintf(stderr,"Missing record of parents!\n");
    exit(1);
  }

  if (!have_trees && logg.data['R']!=0)
  { fprintf(stderr,"Missing record of divergence times!\n");
    exit(1);
  }

  if (!have_trees && (have_latent || have_locations))
  { fprintf(stderr,"Missing records describing tree structures\n");
    exit(1);
  }

  if (m!=0 && (m->type!='R' || m->noise.alpha[2]!=0)
   && !have_latent && !have_locations
   && N_train>0)
  { fprintf(stderr,
"Need latent vectors or node locations to generate cases from non-Gaussian model\n");
    exit(1);
  }

  /* Read training data, if any. */
  
  data_spec = logg.data['D'];

  if (data_spec!=0) 
  { dft_data_read (1, 0, dft, m);
  }
  else
  { N_train = 0;
  }

  if (N_train>0 && !have_trees)
  { fprintf(stderr,"Missing records describing trees for training cases\n");
    exit(1);
  }

  n_tot = N_train+n_cases;

  /* Check sizes of records holding state info. */

  dft_check_sizes(&logg,dft,m,N_train);

  /* Allocate space, and copy over existing data. */

  parents   = chk_alloc(1,dft_parents_size(dft,n_tot));
  divt      = chk_alloc(1,dft_divt_size(dft,n_tot));
  latent    = chk_alloc(1,dft_latent_size(dft,n_tot));
  locations = chk_alloc(1,dft_locations_size(dft,n_tot));

  for (dt = 0; dt<N_trees; dt++)
  { terminals[dt] = chk_alloc (n_tot*N_targets, sizeof (double));
  }

  nodes = chk_alloc (N_trees*n_tot, sizeof *nodes);

  if (N_train>0)
  { 
    divt0 = logg.data['T'];
    j = 0; j0 = 0;
    for (dt = 0; dt<N_trees; dt++)
    { for (i = 0; i<N_train-1; i++)
      { divt[j] = divt0[j0];
        j += 1; j0 += 1;
      }
      j += n_cases;
    }

    parents0 = logg.data['R'];
    j = 0; j0 = 0;
    for (dt = 0; dt<N_trees; dt++)
    { for (i = 0; i<2*N_train; i++)
      { parents[j+n_cases] = parents0[j0];
        j += 1; j0 += 1;
      }
      j += 2*n_cases;
    }

    if (have_latent)
    { latent0 = logg.data['L'];
      for (i = 0; i<N_train; i++)
      { for (t = 0; t<N_targets; t++)
        { latent[i*N_targets+t] = latent0[i*N_targets+t];
        }
      }
    }

    if (have_locations)
    { locations0 = logg.data['N'];
      j = 0; j0 = 0;
      for (dt = 0; dt<N_trees; dt++)
      { for (i = 0; i<N_train-1; i++)
        { for (t = 0; t<N_targets; t++) 
          { locations[j] = locations0[j0];
            j += 1; j0 += 1;
          }
        }
        j += n_cases*N_targets;
      }
    }
  }

  /* Open the output file. */

  output = fopen(output_file,"w");
  if (output==NULL)
  { fprintf(stderr,"Can't create output file (%s)\n",output_file);
    exit(1);
  }

  /* Create nodes and set up state. */

  for (dt = 0; dt<N_trees; dt++)
  { first[dt] = dft_conv_tree (parents+dt*(2*n_tot)+n_tot-1, 
                               nodes+dt*n_tot, N_train);
  }

  dft_setup_state (dft, m, st, hyp, parents, divt, locations, nodes, n_tot);

  /* Generate latent vectors and locations for training cases, if not here. */

  if (!have_locations && N_train>0)
  { if (dft->N_trees>1)
    { fprintf(stderr,
      "Can't sample for absent locations in when there's more than one tree\n");
      exit(1);
    }
    dft_gibbs_locations (dft, m, st, 1, have_latent ? latent : 0);
  }

  if (!have_latent && N_train>0)
  { dft_sample_latent (dft, m, st, 1, latent);
  }

  /* Generate locations for terminal nodes. */

  dft_sample_terminals (dft, st, latent, terminals);

  /* Generate the new cases, writing them to the output file. */

  rand_seed (random_seed);

  for (N = N_train; N<n_tot; N++)
  {
    /* Initialize latent vector for new case to zero. */

    for (t = 0; t<N_targets; t++)
    { latent[N*N_targets+t] = 0; 
    }

    /* Add one more case to each of the trees, incrementing latent vector. */

    for (dt = 0; dt<N_trees; dt++)
    { 
      double *vp, *vc, *vn, *v, *w;
      double s;
      int x;

      v = locations + (n_tot-1)*N_targets*dt;
      w = terminals[dt];

      if (N==0)
      {
        first[dt] = 1;
        st[dt].parents[1] = 0;

        for (t = 0; t<N_targets; t++)
        { w[N*N_targets+t] = st[dt].d_SD[t]*rand_gaussian();
        }
      }
      else
      {
        dft_gen_path (hyp, st, dt, first[dt], &x, &s, (int *) 0);
  
        first[dt] = dft_add_node (st[dt].parents, st[dt].nodes, -N, N+1, x);
        st[dt].divt[N] = s;

        vp = st[dt].parents[-N]==0 ? 0 : v + (-st[dt].parents[-N]-1)*N_targets;
        vc = x>0 ? w+(x-1)*N_targets : v+(-x-1)*N_targets;
        vn = v + (N-1)*N_targets;

        for (t = 0; t<N_targets; t++)
        { double am, bm, av, bv, sd;
          sd = st[dt].d_SD[t];
          am = vp==0 ? 0 : vp[t];
          av = (sd*sd) * (vp==0 ? s : s - st[dt].divt[-st[dt].parents[-N]]);
          bm = vc[t];
          bv = (sd*sd) * (x>0 ? 1-s : st[dt].divt[-x] - s);
          vn[t] = av+bv<=0 ? (am+bm)/2  /* check for this just in case */
                : (am*bv+bm*av)/(av+bv) + rand_gaussian()*sqrt((av*bv)/(av+bv));
          w[N*N_targets+t] = vn[t] + sqrt(1-s)*sd*rand_gaussian();
        }
      }
  
      for (t = 0; t<N_targets; t++)
      { latent[N*N_targets+t] += w[N*N_targets+t];
        if (isnan(latent[N*N_targets+t])) abort();
      }
    }

    /* Generate target values based on the latent vector */

    for (t = 0; t<N_targets; t++)
    { if (m==0)
      { if (!only_tree) fprintf (output, " %+.6e", latent[N*N_targets+t]);
      }
      else if (m->type=='R')
      { double n, x;
        n = prior_pick_sigma (st->noise[t], m->noise.alpha[2]);
        x = latent[N*N_targets+t] + n*rand_gaussian();
        if (!only_tree) fprintf (output, " %+.6e", x);
      }
      else if (m->type=='B')
      { int b;
        b = rand_uniform() < 1/(1+exp(-latent[N*N_targets+t]));
        if (!only_tree) fprintf (output, " %d", b);
      }
      else 
      { abort();
      }
      if (only_tree) 
      { fprintf (output, " %+.6e", terminals[only_tree-1][N*N_targets+t]);
      }
    }

    fprintf(output,"\n");
  }

  /* Print information on standard output, if asked to. */

  if (show_hypers)
  { dft_print_hypers (dft, m, hyp);
  }

  if (show_params)
  { dft_print_params (dft, hyp);
  }

  if (show_latent)
  { dft_print_latent (dft, latent, n_tot, bare);
  }

  if (show_trees)
  { dft_print_trees (dft, st, n_tot, bare);
  }

  if (show_graph)
  { dft_graph_trees (dft, st, n_tot, graph_time, bare);
  }

  if (show_nodes)
  { dft_print_nodes (dft, st, n_tot, include_locations, bare);
  }

  /* Close log file that was read from. */

  log_file_close(&logf);

  /* Write information to new log file, if asked. */

  if (new_log_file)
  {
    logf.file_name = new_log_file;
    log_file_open (&logf, 1);
    log_file_last (&logf);

    index = logf.at_end ? 0 : logf.header.index+1;
    
    logf.header.type = 'S';
    logf.header.index = index;
    logf.header.size = dft_hyper_size(dft,m);
    log_file_append (&logf, hyp);

    logf.header.type = 'T';
    logf.header.index = index;
    logf.header.size = dft_divt_size(dft,n_tot);
    log_file_append (&logf, divt);

    logf.header.type = 'R';
    logf.header.index = index;
    logf.header.size = dft_parents_size(dft,n_tot);
    log_file_append (&logf, parents);
  
    logf.header.type = 'L';
    logf.header.index = index;
    logf.header.size = dft_latent_size(dft,n_tot);
    log_file_append (&logf, latent);
  
    logf.header.type = 'N';
    logf.header.index = index;
    logf.header.size = dft_locations_size(dft,n_tot);
    log_file_append (&logf, locations);

    log_file_close(&logf);
  }

  exit(0);
}
