/* NET-DVAR.C - Program to find variance of function difference. */

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
  net_arch   *a;
  net_flags  *flgs;
  net_priors *p;

  net_params params, *w = &params;
  net_values values1, *v1 = &values1;
  net_values values2, *v2 = &values2;

  double f1, f2;
  double *dvar;

  int hidden;
  int N_points;

  log_file logf;
  log_gobbled logg;

  int lindex, hindex, index_mod;
  double centre, low_width, high_width;
  int n_pairs;

  double width;
  int i, j;

  /* Look at arguments. */

  hidden = 0;
  if (argc>0 && strcmp(argv[argc-1],"hidden")==0)
  { hidden = 1;
    argc -= 1;
  }

  if (argc!=7 && argc!=8 || strcmp(argv[3],"/")!=0)
  { usage();
  }

  centre = atof(argv[4]);
  n_pairs = atoi(argv[5]);
  low_width = argc==7 ? 0 : atof(argv[6]);
  high_width = atof(argv[argc-1]);

  if (argc==8 && low_width<=0 || high_width<=0 || n_pairs<=0) usage();

  parse_range (argv[2], &lindex, &hindex, &index_mod);

  if (hindex<0)
  { hindex = 1000000000;
  }

  /* Open log file and read network architecture and priors. */

  logf.file_name = argv[1];

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  net_record_sizes(&logg);

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }
  
  if ((a = logg.data['A'])==0)
  { fprintf(stderr,"No architecture specification in log file\n");
    exit(1);
  }

  flgs = logg.data['F'];

  w->total_params = net_setup_param_count(a,flgs);

  logg.req_size['W'] = w->total_params * sizeof(net_param);

  if ((p = logg.data['P'])==0)
  { fprintf(stderr,"No prior specification in log file\n");
    exit(1);
  }

  if (hidden && a->N_layers==0)
  { fprintf(stderr,"Network doesn't have hidden units\n");
    exit(1);
  }

  /* Allocate space for network values and variances. */

  net_setup_value_pointers
    (v1, chk_alloc (net_setup_value_count(a), sizeof(net_value)), a);
  net_setup_value_pointers
    (v2, chk_alloc (net_setup_value_count(a), sizeof(net_value)), a);

  dvar = chk_alloc (n_pairs+1, sizeof *dvar);

  /* Go through all the networks, collecting statistics. */

  for (i = 0; i<=n_pairs; i++) dvar[i] = 0;
  for (i = 1; i<a->N_inputs; i++) v1->i[i] = v2->i[i] = 0;
  N_points = 0;

  for (;;)
  {
    /* Skip to next desired index, or to end of range. */

    while (!logf.at_end && logf.header.index<=hindex
     && (logf.header.index<lindex || (logf.header.index-lindex)%index_mod!=0))
    { log_file_forward(&logf);
    }

    if (logf.at_end || logf.header.index>hindex)
    { break;
    }

    /* Gobble up records for this index. */

    log_gobble(&logf,&logg);
   
    /* Use the network with this index, if such a network exists. */

    if (logg.index['W']==logg.last_index)
    {
      w->param_block = logg.data['W'];

      net_setup_param_pointers (w, a, flgs);

      if (hidden)
      { 
        for (i = 0; i<=n_pairs; i++)
        { width = low_width==0 ? high_width*i/n_pairs
                : low_width*pow(high_width/low_width,(double)i/n_pairs);
          v1->i[0] = centre + width/2;
          net_func (v1, 0, a, flgs, w);
          v2->i[0] = centre - width/2;
          net_func (v2, 0, a, flgs, w);
          for (j = 0; j<a->N_hidden[a->N_layers-1]; j++)
          { f1 = v1->h[a->N_layers-1][j];
            f2 = v2->h[a->N_layers-1][j];
            dvar[i] += (f1-f2)*(f1-f2);
          }
        }

        N_points += a->N_hidden[a->N_layers-1];

      }
      else
      {
        for (i = 0; i<=n_pairs; i++)
        { width = low_width==0 ? high_width*i/n_pairs
                : low_width*pow(high_width/low_width,(double)i/n_pairs);
          v1->i[0] = centre + width/2;
          net_func (v1, 0, a, flgs, w);
          v2->i[0] = centre - width/2;
          net_func (v2, 0, a, flgs, w);
          f1 = v1->o[0];
          f2 = v2->o[0];
          dvar[i] += (f1-f2)*(f1-f2);
        }

        N_points += 1;
      }
    }
  }

  if (N_points==0) 
  { fprintf(stderr,"No networks found in specified range\n");
    exit(1);
  }

  /* Print variances. */

  for (i = 0; i<=n_pairs; i++)
  { width = low_width==0 ? high_width*i/n_pairs
          : low_width*pow(high_width/low_width,(double)i/n_pairs);
    printf ("%20.6le %20.6le\n", width, dvar[i]/N_points);
  }
  
  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf(stderr,
   "Usage: net-dvar log-file range / centre n-pairs [ low ] high [ \"hidden\" ]\n");
  exit(1);
}
