/* MOL-QUANTITIES.C - Module defining molecular dynamics quantities. */

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
#include "quantities.h"
#include "mc.h"
#include "mol.h"


static mol_spec *ms;


/* INITIALIZE AFTER FIRST RECORDS READ. */

void mol_initialize
( log_gobbled *logg
)
{ 
  ms = logg->data['M'];

  if (ms==0)
  { fprintf(stderr,"Molecular dynamics records not present\n");
    exit(1);
  }

  if (logg->actual_size['C']!=((ms->len_pres<0)+ms->D*ms->N)*sizeof(mc_value))
  { fprintf(stderr,"Bad size for state records\n");
    exit(1);
  }
}


/* INDICATE WHAT QUANTITIES ARE AVAILABLE FROM THIS MODULE. */

void mol_available 
( quantities_described qd,
  log_gobbled *logg
)
{ 
  char letter;
  int mod;
  int v;

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    mod = qd[v].modifier;

    if (letter && qd[v].available==0)
    {
      if ((letter=='x' || letter=='y' || letter=='z') 
            && (mod==-1 || mod==1) && qd[v].low!=-1)
      { qd[v].available = letter=='y' && ms->D<2 || letter=='z' && ms->D<3
                           || qd[v].high>=ms->N ? -1 : 1;
        if (qd[v].high==-1) qd[v].high = ms->N-1;
      }
      else if (strchr("PpVvUuWwOo",letter) && mod==-1 
            && qd[v].low==-1 && qd[v].high==-1)
      { qd[v].available = 1;
      }
      else if ((letter=='n' || letter=='N') && mod==-1)
      { qd[v].available = qd[v].high<ms->N && qd[v].low<ms->N ? 1 : -1;
        if (qd[v].high==-1 && qd[v].low!=-1) qd[v].high = ms->N-1;
      }
      else if (letter=='d' && mod!=-1 && (qd[v].low!=-1 || qd[v].high!=-1))
      { qd[v].available = mod<ms->N && qd[v].high<ms->N && qd[v].low<ms->N 
                           ? 1 : -1;
        if (qd[v].high==-1) qd[v].high = ms->N-1;
        if (qd[v].low==-1)  qd[v].low = 0;
      }
    }
  }
}


/* EVALUATE QUANTITIES KNOWN TO THIS MODULE. */

void mol_evaluate 
( quantities_described qd, 
  quantities_held *qh,
  log_gobbled *logg
)
{ 
  double *coords;
  double len;

  int mod, low, high;
  char letter;
  int v, i, j, k, l, h;

  if (logg->data['C']==0 || logg->index['C']!=logg->last_index)
  { fprintf(stderr,"State record missing\n"); 
    return;
  }
  
  coords = logg->data['C'];
  len = dlen(ms,coords);

  for (v = 0; v<Max_quantities; v++)
  {
    letter = qd[v].letter;
    mod = qd[v].modifier;
    low = qd[v].low;
    high = qd[v].high;

    if (letter && !qh->updated[v])
    {
      if (letter=='x' || letter=='y' || letter=='z')
      { 
        for (i = low; i<=high; i++)
        { if (mod==1)
          { qh->value[v][i-low] = len * wrap(0.5+coords[ms->D*i+(letter-'x')]);
          }
          else
          { qh->value[v][i-low] = len * wrap (coords[ms->D*i+(letter-'x')]);
          }
        }

        qh->updated[v] = 1;
      }

      else if (mod==-1 &&
         (letter=='V' || letter=='v' || letter=='P' || letter=='p'))
      {
        double gri[3], grj[3], c;

        *qh->value[v] = 0;

        c = ms->width * pow (2.0, 1/6.0);  /* Point of minimum potential */

        for (i = 0; i<ms->N; i++)
        { for (j = i+1; j<ms->N; j++)
          { double d, f;
            d = sqrt (squared_distance (ms->D, len, coords, i, j));
            (void) energy_and_gradient (coords, 1, ms, i, j, gri, 0, 0);
            f = 0;
            for (k = 0; k<ms->D; k++)
            { f += (gri[k] * gri[k]) / (len * len);
            }
            f = sqrt(f);
            if (d>c) f = -f;
            *qh->value[v] += f * d;
          }
        }

        *qh->value[v] /= ms->D;

        if (letter=='P' || letter=='p')
        { *qh->value[v] += ms->N;
          *qh->value[v] /= pow(len,(double)ms->D);
        }

        if (letter=='p' || letter=='v')
        { *qh->value[v] /= ms->scale;
        }

        if (letter=='p')
        { *qh->value[v] *= pow(ms->width,ms->D);
        }

        qh->updated[v] = 1;
      }

      else if (letter=='w')
      { *qh->value[v] = len;
        qh->updated[v] = 1;
      }
      else if (letter=='W')
      { *qh->value[v] = pow(len,ms->D);
        qh->updated[v] = 1;
      }
      else if (letter=='O')
      { *qh->value[v] = ms->N / pow(len,ms->D);
        qh->updated[v] = 1;
      }
      else if (letter=='o')
      { *qh->value[v] = pow(ms->width,ms->D) * ms->N / pow(len,ms->D);
        qh->updated[v] = 1;
      }

      else if (letter=='U' || letter=='u')
      { *qh->value[v] = 0;
        for (i = 0; i<ms->N; i++)
        { for (j = i+1; j<ms->N; j++)
          { *qh->value[v] += energy_and_gradient (coords, 1, ms, i, j, 0, 0, 0);
          }
        }
        *qh->value[v] /= ms->scale;
        if (letter=='u') 
        { *qh->value[v] /= ms->N;
        }
        qh->updated[v] = 1;
      }

      else if (letter=='n' || letter=='N')
      { 
        double minmax, minmaxa, d2;
        l = low==-1 ? 0 : low;
        h = high==-1 ? ms->N-1 : high;
        minmaxa = letter=='n' ? ms->D*len*len : 0;
        for (i = l; i<=h; i++)
        { minmax = letter=='n' ? ms->D*len*len : 0;
          for (j = 0; j<ms->N; j++)
          { if (j!=i)
            { d2 = squared_distance(ms->D,len,coords,i,j);
              if (letter=='n' && d2<minmax || letter=='N' && d2>minmax) 
              { minmax = d2;
              }
            }
          }
          if (low==-1)
          { if (letter=='n' && minmax<minmaxa || letter=='N' && minmax>minmaxa) 
            { minmaxa = minmax;
            }
          }
          else 
          { qh->value[v][i-low] = sqrt(minmax);
          }
        }
        if (low==-1)
        { *qh->value[v] = sqrt(minmaxa);
        }
        qh->updated[v] = 1;
      }

      else if (mod!=-1 && letter=='d')
      {
        for (i = low; i<=high; i++)
        { qh->value[v][i-low] = 
            sqrt(squared_distance(ms->D,len,coords,mod,i));
        }
        qh->updated[v] = 1;
      }
    }
  }
}
