/* MOL-DISPLAY.C - Program to display state of molecular dynamics simulation. */

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


/* MAIN PROGRAM. */

main
( int argc,
  char **argv
)
{
  mol_spec *ms;
  double *coords;

  log_file logf;
  log_gobbled logg;

  int index, shift, plen;
  int j, i, iD;
  double len;

  /* Look at arguments. */

  index = -1;
  shift = 0;
  plen = 0;

  if (argc>1 && strcmp(argv[1],"-s")==0)
  { shift = 1;
    argc -= 1;
    argv += 1;
  }
  else if (argc>1 && strcmp(argv[1],"-l")==0)
  { plen = 1;
    argc -= 1;
    argv += 1;
  }

  if (argc!=2 && argc!=3 
   || argc>2 && (index = atoi(argv[2]))<=0 && strcmp(argv[2],"0")!=0) 
  { fprintf (stderr, "Usage: mol-display [ -s | -l ] log-file [ index ]\n");
    exit(1);
  }

  logf.file_name = argv[1];

  /* Open log file and read specification. */

  log_file_open (&logf, 0);

  log_gobble_init(&logg,0);
  logg.req_size['M'] = sizeof *ms;

  while (!logf.at_end && logf.header.index<0)
  { log_gobble(&logf,&logg);
  }

  ms = logg.data['M'];
  
  if (ms==0)
  { fprintf(stderr,"No molecular dynamics specification in log file\n");
    exit(1);
  }

  if (plen && ms->len_pres>0)
  { fprintf(stderr,
 "The -l option is allowed only when using the NPT ensemble (variable length)\n"
    );
    exit(1);
  }

  /* Read the desired state from the log file. */

  if (index<0)
  { 
    log_gobble_last(&logf,&logg);

    if (logg.last_index<0)
    { fprintf(stderr,"No molecular dynamics state in log file\n");
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

  if (logg.index['C']!=index)
  { fprintf(stderr,"No molecular dynamics state stored with that index\n");
    exit(1);
  }

  if (logg.actual_size['C']!=(ms->D*ms->N+(ms->len_pres<0))*sizeof(mc_value))
  { fprintf(stderr,"Bad size for state records\n");
    exit(1);
  }

  coords = logg.data['C'];

  /* Print length of dimensions if -l specified. */

  if (plen)
  { printf(" %10.6f\n",exp(coords[ms->D*ms->N]));
  }

  /* Otherwise, print the coordinates of the particles. */
  
  if (!plen)
  {
    len = dlen(ms,coords);

    for (i = 0; i<ms->N; i++)
    { 
      iD = i * ms->D;
  
      for (j = 0; j<ms->D; j++)
      { printf (" %10.6f", len * wrap(0.5*shift + coords[iD+j]));
      }

      printf("\n");
    }
  }

  exit(0);
}
