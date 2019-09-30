/* DATA-SPEC.C - Program for specifying data sets for training and testing. */

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
#include "numin.h"
#include "data.h"
#include "log.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  /* The following variables are static so that the total size of local 
     variables on the stack is kept down, just in case a compiler might have 
     problems with more that 64K bytes in a stack frame. */

  static numin_source ns;

  static double iv[Max_inputs];
  static double tv[Max_targets];

  static double imean[Max_inputs];
  static double ivar[Max_inputs];

  static double tmean[Max_targets];
  static double tvar[Max_targets];
  static int tn[Max_targets];

  data_specifications *ds;

  int do_not_read_data;
  int N_inputs, N_targets;

  log_file logf;
  log_gobbled logg;
  int anytrans;
  int i, j, n;
  char **ap;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing data specifications. */

  if (argc==2)
  {
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }

    if (logg.data['D']==0)
    { fprintf(stderr,"No data specifications in log file\n");
      exit(1);
    }

    ds = logg.data['D'];

    N_inputs = ds->N_inputs;
    N_targets = ds->N_targets;

    if (logg.actual_size['D'] != data_spec_size(N_inputs,N_targets))
    { fprintf(stderr,"Data specification records is the wrong size!\n");
      exit(1);
    }

    printf("\nData specifications:\n\n");

    printf("  Number of inputs:  %d\n",N_inputs);
    printf("  Number of targets: %d\n",N_targets);
    printf("\n");

    if (ds->int_target) 
    { printf("  Targets are integers in the range [0,%d)\n",ds->int_target);
      printf("\n");
    }

    anytrans = 0;

    for (i = 0; i<N_inputs; i++)
    { if (ds->trans[i].take_log 
       || ds->trans[i].data_shift
       || ds->trans[i].data_scale
       || ds->trans[i].shift!=0
       || ds->trans[i].scale!=1)
      { anytrans = 1;
        break;
      }
    }

    if (anytrans)
    { 
      printf("  Transformations for input variables:\n\n");

      for (i = 0; i<N_inputs; i++)
      { printf("    %2d: %s\n",i+1,data_trans_build(ds->trans[i]));
      }
  
      printf("\n");
    }

    anytrans = 0;

    for (i = 0; i<N_targets; i++)
    { if (ds->trans[N_inputs+i].take_log 
       || ds->trans[N_inputs+i].data_shift
       || ds->trans[N_inputs+i].data_scale
       || ds->trans[N_inputs+i].shift!=0
       || ds->trans[N_inputs+i].scale!=1)
      { anytrans = 1;
        break;
      }
    }

    if (anytrans)
    { 
      printf("  Transformations for target variables:\n\n");

      for (i = 0; i<N_targets; i++)
      { printf("    %2d: %s\n",i+1,data_trans_build(ds->trans[N_inputs+i]));
      }
  
      printf("\n");
    }

    printf ("  Training inputs:  %s\n", ds->train_inputs);
    printf ("  Training targets: %s\n", ds->train_targets);
    printf("\n");

    if (ds->test_inputs[0]!=0)
    { printf ("  Test inputs:  %s\n", ds->test_inputs);
      if (ds->test_targets[0]!=0)
      { printf ("  Test targets: %s\n", ds->test_targets);
      }
      printf("\n");
    }
    
    log_file_close(&logf);
  
    exit(0);
  }

  /* Otherwise, look at remaining arguments, up to transformations. */

  ap = argv+2;

  if (*ap==0 || (N_inputs = atoi(*ap++))<0 
   || N_inputs==0 && strcmp(*(ap-1),"0")!=0) usage();
  if (*ap==0 || (N_targets = atoi(*ap++))<0
   || N_targets==0 && strcmp(*(ap-1),"0")!=0) usage();

  if (N_inputs>Max_inputs)
  { fprintf(stderr,"Too many input values (max %d)\n",Max_inputs);
    exit(1);
  }
  if (N_targets>Max_targets)
  { fprintf(stderr,"Too many target values (max %d)\n",Max_targets);
    exit(1);
  }

  ds = chk_alloc (1, data_spec_size(N_inputs,N_targets));

  ds->N_inputs = N_inputs;
  ds->N_targets = N_targets;

  ds->int_target = 0;
  ds->test_inputs[0] = 0;
  ds->test_targets[0] = 0;

  if (*ap!=0 && strcmp(*ap,"/")!=0 && strcmp(*ap,"-n")!=0)
  { if ((ds->int_target = atoi(*ap++))<=0) usage();
  }

  do_not_read_data = 0;
  if (*ap!=0 && strcmp(*ap,"-n")==0)
  { do_not_read_data = 1;
    ap += 1;
  }

  if (*ap==0 || strcmp(*ap++,"/")!=0) usage();

  if (*ap==0) usage();
  strcpy(ds->train_inputs,*ap++);

  if (*ap==0) usage();
  strcpy(ds->train_targets,*ap++);

  if (*ap!=0 && strcmp(*ap,"/")!=0)
  { strcpy(ds->test_inputs,*ap++);
    if (*ap!=0 && strcmp(*ap,"/")!=0)
    { strcpy(ds->test_targets,*ap++);
    }
  }

  /* Set transformations to the identity. */

  for (i = 0; i<N_inputs; i++)
  { ds->trans[i] = data_trans_parse("I");
  }

  for (i = 0; i<N_targets; i++)
  { ds->trans[ds->N_inputs+i] = data_trans_parse("I");
  }

  /* Look at transformation arguments (if any). */

  if (*ap!=0 && strcmp(*ap,"/")==0)
  { 
    ap += 1;

    for (i = 0; 
         i<ds->N_inputs && *ap!=0 && strcmp(*ap,"/")!=0 && strcmp(*ap,"...")!=0;
         i++)
    { ds->trans[i] = data_trans_parse(*ap++);
    }

    if (*ap!=0 && strcmp(*ap,"...")==0)
    { if (i==0) usage();
      for ( ; i<ds->N_inputs; i++)
      { ds->trans[i] = ds->trans[i-1];
      }
      ap += 1;
    }

    if (*ap!=0 && strcmp(*ap,"/")!=0)
    { fprintf(stderr,"Too many input variable transformations\n");
      exit(1);
    }

    if (*ap!=0 && strcmp(*ap,"/")==0)
    { 
      ap += 1;

      for (i = 0; 
        i<ds->N_targets && *ap!=0 && strcmp(*ap,"/")!=0 && strcmp(*ap,"...")!=0;
        i++)
      { ds->trans[ds->N_inputs+i] = data_trans_parse(*ap++);
        if (ds->int_target && (ds->trans[ds->N_inputs+i].data_shift 
                                || ds->trans[ds->N_inputs+i].data_scale))
        { fprintf(stderr,
          "Can't use data-derived target transformation for integer targets\n");
          exit(1);
        }
      }

      if (*ap!=0 && strcmp(*ap,"...")==0)
      { if (i==0) usage();
        for ( ; i<ds->N_targets; i++)
        { ds->trans[ds->N_inputs+i] = ds->trans[ds->N_inputs+i-1];
        }
        ap += 1;
      }

      if (*ap!=0 && strcmp(*ap,"/")!=0)
      { fprintf(stderr,"Too many target transformations\n");
        exit(1);
      }
    }
  }

  if (*ap!=0) usage();

  /* Read training data, and compute any data-dependent shift and scale 
     amounts. */

  if (!do_not_read_data)
  {
    numin_spec (&ns, "data@1,0", 1);
    numin_spec (&ns, ds->train_inputs, ds->N_inputs);
  
    n = numin_start(&ns);
  
    for (j = 0; j<ds->N_inputs; j++)
    { imean[j] = 0;
      ivar[j] = 0;
    }
  
    for (i = 0; i<n; i++)
    { numin_read(&ns,iv);
      for (j = 0; j<ds->N_inputs; j++)
      { if (isnan(iv[j])) 
        { fprintf(stderr,"Missing values not allowed for inputs\n");
          exit(1);
        }
        iv[j] = data_trans (iv[j], ds->trans[j]);
        imean[j] += iv[j];
        ivar[j] += iv[j]*iv[j];
      }
    }
  
    for (j = 0; j<ds->N_inputs; j++)
    { imean[j] /= n;
      ivar[j] /= n;
      ivar[j] -= imean[j]*imean[j];
      if (ds->trans[j].data_shift)
      { ds->trans[j].shift = -imean[j];
      }
      if (ds->trans[j].data_scale)
      { ds->trans[j].scale = ivar[j]<=0 ? 1 : 1/sqrt(ivar[j]);
      }
    }
  
    numin_close(&ns);
  
    numin_spec (&ns, ds->train_targets, ds->N_targets);
  
    if (numin_start(&ns)!=n)
    { fprintf(stderr,
        "Number of training input cases doesn't match number of targets\n");
      exit(1);
    }
  
    for (j = 0; j<ds->N_targets; j++)
    { tmean[j] = 0;
      tvar[j] = 0;
      tn[j] = 0;
    }
  
    for (i = 0; i<n; i++)
    { numin_read(&ns,tv);
      for (j = 0; j<ds->N_targets; j++)
      { if (isnan(tv[j])) continue;
        tv[j] = data_trans (tv[j], ds->trans[ds->N_inputs+j]);
        tmean[j] += tv[j];
        tvar[j] += tv[j]*tv[j];
        tn[j] += 1;
        if (ds->int_target)
        { if (tv[j]!=(int)tv[j] || tv[j]<0 || tv[j]>=ds->int_target)
          { fprintf(stderr,"Training target out of bounds or not integer: %f\n",
                    tv[j]);
            exit(1);
          }
        }
      }
    }
  
    for (j = 0; j<ds->N_targets; j++)
    { tmean[j] /= tn[j];
      tvar[j] /= tn[j];
      tvar[j] -= tmean[j]*tmean[j];
      if (ds->trans[ds->N_inputs+j].data_shift)
      { ds->trans[ds->N_inputs+j].shift = -tmean[j];
      }
      if (ds->trans[ds->N_inputs+j].data_scale)
      { ds->trans[ds->N_inputs+j].scale = tvar[j]<=0 ? 1 : 1/sqrt(tvar[j]);
      }
    }
  
    numin_close(&ns);
  
    fprintf(stderr,"Number of training cases: %d\n",n);
  
    /* Read test data. */
  
    if (*ds->test_inputs)
    {
      numin_spec (&ns, "data@1,0", 1);
      numin_spec (&ns, ds->test_inputs, ds->N_inputs);
    
      n = numin_start(&ns);
    
      for (i = 0; i<n; i++)
      { numin_read(&ns,iv);
        for (j = 0; j<ds->N_targets; j++)
        { if (isnan(iv[j])) 
          { fprintf(stderr,"Missing values not allowed for inputs\n");
            exit(1);
          }
          iv[j] = data_trans (iv[j], ds->trans[j]);
        } /* Above done just so errors will be reported at this time */
      }
  
      numin_close(&ns);
   
      if (*ds->test_targets)
      {  
        numin_spec (&ns, ds->test_targets, ds->N_targets);
      
        if (numin_start(&ns)!=n)
        { fprintf(stderr,
            "Number of test input cases doesn't match number of targets\n");
        }
      
        for (i = 0; i<n; i++)
        { numin_read(&ns,tv);
          for (j = 0; j<ds->N_targets; j++)
          { if (isnan(tv[j])) continue;
            tv[j] = data_trans (tv[j], ds->trans[ds->N_inputs+j]);
          }
          if (ds->int_target)
          { for (j = 0; j<ds->N_targets; j++)
            { if (isnan(tv[j])) continue;
              if (tv[j]!=(int)tv[j] || tv[j]<0 || tv[j]>=ds->int_target)
              { fprintf(stderr,"Test target out of bounds or not integer: %f\n",
                        tv[j]);
                exit(1);
              }
            }
          }
        }
  
        numin_close(&ns);
      }
  
      fprintf(stderr,"Number of test cases: %d\n",n);
    }
  }

  /* Open log file. */

  log_file_open (&logf, 1);

  /* Write data specifications to log file. */

  logf.header.type = 'D';
  logf.header.index = -1;
  logf.header.size = data_spec_size(N_inputs,N_targets);
  log_file_append (&logf, ds);

  log_file_close (&logf);

  exit(0);
}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{ 
  fprintf(stderr,
   "Usage: data-spec log-file N-inputs N-targets [ int-target ] [ -n ]\n");
  fprintf(stderr,
   "                   / train-inputs train-targets [ test-inputs [ test-targets ] ]\n");
  fprintf(stderr,
   "                 [ / { input-trans } [ / { target-trans } ] ]\n");
  fprintf(stderr,
   "   or: data-spec log-file (to display stored data specifications)\n");

  exit(1);
}
