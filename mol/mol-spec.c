/* MOL-SPEC.C - Specify molecular dynamics system. */

/* Copyright (c) 1995-2003 by Radford M. Neal 
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
#include "mc.h"
#include "mol.h"


static void usage(void);


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
) 
{
  static mol_spec mspec; /* static so unused fields are intialized to zero */
  mol_spec *ms = &mspec;

  log_file logf;
  log_gobbled logg;

  char **ap, *ep;
  char junk;

  /* Look for log file name. */

  if (argc<2) usage();

  logf.file_name = argv[1];

  /* See if we are to display existing specifications. */

  if (argc==2)
  {
    /* Open log file and gobble up initial records. */
  
    log_file_open(&logf,0);

    log_gobble_init(&logg,0);
    logg.req_size['M'] = sizeof *ms;

    if (!logf.at_end && logf.header.index==-1)
    { log_gobble(&logf,&logg);
    }

    /* Display specification. */  

    if ((ms = logg.data['M'])==0)
    { fprintf(stderr,"No molecular dynamics specification in log file\n");
      exit(1);
    }

    printf("\n");
    printf("%d-D MOLECULAR DYNAMICS MODEL, %s ENSEMBLE\n\n",
      ms->D, ms->len_pres>0 ? "NVT" : "NPT");

    printf("Scale parameter for energy: %.4f\n",ms->scale);
    printf("Width parameter for energy: %.4f\n",ms->width);

    printf("\n");
    if (ms->max_pair_energy>0) 
    { printf("Maximum pair energy: %.1f\n",ms->max_pair_energy);
      printf("\n");
    }

    printf("Number of molecules:        %d\n",ms->N);
    if (ms->len_pres>0)
    { printf("Length in each dimension:   %.4f\n",ms->len_pres);
      printf("Actual density:             %.4f\n", 
                 ms->N/pow(ms->len_pres,(double)ms->D));
      printf("Reduced density:            %.4f\n", 
                 ms->N/pow(ms->len_pres/ms->width,(double)ms->D));
    }
    else
    { printf("Actual pressure:            %.4f\n",
                 -ms->len_pres);
      printf("Reduced pressure:           %.4f\n",
                 -ms->len_pres*pow(ms->width,(double)ms->D)/ms->scale);
    }

    printf("Reduced temperature:        %.4f\n", 1.0/ms->scale);
    printf("\n");

    if (ms->inv_temp!=0)
    { printf(
         "Inverse temperature to initialize volume for AIS, HIS, MIS: %.6lf\n",
          ms->inv_temp);
      printf("\n");
    }

    log_file_close(&logf);

    exit(0);
  }

  /* Otherwise, look at arguments. */

  ap = argv+2;

  if (*ap==0 || sscanf(*ap++,"%d%c",&ms->D,&junk)!=1 || ms->D<1 || ms->D>3
   || *ap==0 || sscanf(*ap++,"%lf%c",&ms->scale,&junk)!=1 || ms->scale<0
   || *ap==0 || sscanf(*ap++,"%lf%c",&ms->width,&junk)!=1 || ms->width<=0)
  { usage();
  }

  if (*ap!=0 && strcmp(*ap,"NVT")!=0 && strcmp(*ap,"NPT")!=0)
  { if (sscanf(*ap++,"%lf%c",&ms->max_pair_energy,&junk)!=1)
    { usage();
    }
  }

  ep = *ap++;
  if (ep==0 || strcmp(ep,"NVT")!=0 && strcmp(ep,"NPT")!=0)
  { usage();
  }

  if (*ap==0 || sscanf(*ap++,"%d%c",&ms->N,&junk)!=1 || ms->N<1
   || *ap==0 || sscanf(*ap++,"%lf%c",&ms->len_pres,&junk)!=1 || ms->len_pres<=0)
  { usage();
  }

  if (strcmp(ep,"NPT")==0)
  { ms->len_pres = -ms->len_pres;
    if (*ap!=0)
    { if (sscanf ((**ap=='/')+*ap, "%lf%c", &ms->inv_temp,&junk) != 1
       || ms->inv_temp<=0)
      { usage();
      }
      if (**ap=='/') 
      { ms->inv_temp = 1 / ms->inv_temp;
      }      
      ap += 1;
    }
  }

  if (*ap!=0) 
  { usage();
  }

  /* Create log file and write records. */

  log_file_create(&logf);

  logf.header.type = 'M';
  logf.header.index = -1;
  logf.header.size = sizeof *ms;
  log_file_append(&logf,ms);

  log_file_close(&logf);

  exit(0);

}


/* DISPLAY USAGE MESSAGE AND EXIT. */

static void usage(void)
{
  fprintf(stderr,
"Usage: mol-spec log-file dim scale width [ maxpe ] \"NVT\" N length\n");
  fprintf(stderr,
"   or: mol-spec log-file dim scale width [ maxpe ] \"NPT\" N pressure [ inv-temp ]\n");
  fprintf(stderr,
"   or: mol-spec log-file  (to display stored specifications)\n");

  exit(1);
}
