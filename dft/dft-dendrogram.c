/* DFT-DENDROGRAM.C - Program to plot dendrogram for diffusion tree. */

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


/* PRINT USAGE MESSAGE. */

static void usage(void)
{ fprintf (stderr, 
    "Usage: dft-dendrogram [ -r ] log-file index [ tree ] [ label-file ]\n");
  exit(1);
}


/* PRODUCE DENDROGRAM.  Called recursively to do the parts. */

static int dendrogram
( dft_spec *dft,
  dft_state st,
  int dt,
  int x,
  int hoff,
  int voff,
  int max_depth,
  char **labels
)
{ 
  int i, n0, n1, na, nb, ca, cb, va, vb;

  if (x>0)
  { printf("newpath %d %d moveto %d %d lineto stroke %d %d moveto",
            hoff, voff, 50+max_depth*10, voff, 55+max_depth*10, voff-3);
    if (labels)
    { printf(" (%s) show\n",labels[x-1]);
    }
    else
    { printf(" (%d) show\n",x);
    }
    return voff;
  }
  else
  { 
    n0 = npts (st[dt].nodes[-x], 0);
    n1 = npts (st[dt].nodes[-x], 1);
    ca = chld (st[dt].nodes[-x], 0);
    cb = chld (st[dt].nodes[-x], 1);
    na = n0;
    nb = n1;
    if (n0<n1 || n0==n1 && ca>cb)
    { ca = chld (st[dt].nodes[-x], 1);
      cb = chld (st[dt].nodes[-x], 0);
      na = n1;
      nb = n0;
    }

    va = dendrogram (dft, st, dt, ca, hoff+10, voff+10*nb, max_depth, labels);
    vb = dendrogram (dft, st, dt, cb, hoff+10, voff,       max_depth, labels); 

    printf(
      "newpath %d %d moveto %d %d lineto %d %d moveto %d %d lineto stroke\n",
       hoff, voff+10*nb-5, hoff+10, voff+10*nb-5, hoff+10, va, hoff+10, vb);

    return voff+10*nb-5;
  }
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

  int index;
  int dt;
  char *label_file;

  int *parents;	
  double *divt;	
  double *latent;
  double *locations;

  int N_train;

  dft_state st;
  dft_tree_node *nodes;
 
  char junk;
  int max_depth, root, max_label;
  int i;

  int rotate;
  int bblx, bbly, bbux, bbuy;

  FILE *lf;
  char **labels;

  /* Look at arguments. */

  rotate = 0;
  if (argc>1 && strcmp(argv[1],"-r")==0)
  { rotate = 1;
    argc -= 1;
    argv += 1;
  }

  if (argc<3 && argc>5) usage();

  logf.file_name = argv[1];

  if (sscanf(argv[2],"%d%c",&index,&junk)!=1 || index<0)
  { usage();
  }

  label_file = argc>3 ? argv[argc-1] : 0;

  dt = 1;
  if (argc>3)
  { if (sscanf(argv[3],"%d%c",&i,&junk)==1)
    { dt = i;
      if (argc==4)
      { label_file = 0;
      }
    }
  }
  dt -= 1;

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

  if (dt>=dft->N_trees)
  { fprintf(stderr,"There is no tree %d for this model\n",dt+1);
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

    (void) dft_conv_tree (parents+dt*(2*N_train)+N_train-1, 
                          nodes+dt*N_train, N_train);
  }

  dft_setup_state (dft, m, st, h, parents, divt, locations, nodes, N_train);

  /* Read labels. */

  labels = 0;
  max_label = 1 + (int) (0.0001 + log(N_train)/log(10.0));

  if (label_file)
  { lf = fopen(label_file,"r");
    if (lf==NULL)
    { fprintf(stderr,"Can't open file of labels (%s)\n",label_file);
      exit(1);
    }
    labels = chk_alloc (N_train, sizeof *labels);
    for (i = 0; i<N_train; i++)
    { char s[1000];
      if (fscanf(lf,"%999s",s)!=1)
      { fprintf(stderr,"Error reading label file\n");
        exit(1);
      }
      if (strlen(s)==999)
      { fprintf(stderr,"Label is too long\n");
        exit(1);
      }
      labels[i] = chk_alloc(strlen(s)+1,1);
      if (strlen(s)>max_label) max_label = strlen(s);
      strcpy(labels[i],s);
    }
  }

  /* Produce dendrogram. */

  for (root = 1; st[dt].parents[root]!=0; root = st[dt].parents[root]) ;

  max_depth = dft_max_depth(st[dt].parents,N_train);

  printf("%%!PS-Adobe EPSF-3.0\n");
  printf("%%%%Pages: 1\n");

  if (rotate)
  { bblx = 50;
    bbly = 50;
    bbux = 50+N_train*10;
    bbuy = 50+max_depth*10+2+7*max_label;
  }
  else
  { bblx = 50;
    bbly = 50;
    bbux = 50+max_depth*10+2+7*max_label;
    bbuy = 50+N_train*10;
  }

  printf("%%%%BoundingBox: %d %d %d %d\n",bblx,bbly,bbux,bbuy);

  if (0) /* For debugging the bounding box */
  {
    printf(
    "newpath %d %d moveto %d %d lineto %d %d lineto %d %d lineto closepath\n",
     bblx, bbly, bblx, bbuy, bbux, bbuy, bbux, bbly);
    printf("stroke\n");
  }

  printf("%%%%EndComments\n");

  printf("/Courier-Bold findfont 9 scalefont setfont\n");

  if (rotate)
  { printf("%d %d translate -90 rotate %d %d translate\n",
      bblx, bbly, -bbuy, -bbly);
  }

  dendrogram(dft,st,dt,root,50,55,max_depth,labels);

  printf("showpage\n");

  exit(0);
}
