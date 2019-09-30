/* FIND-MIN.C - Find entry with minimum value. */

#include <stdio.h>
#include <stdlib.h>

main()
{
  double t, x, mt, mx;

  if (scanf("%lf%lf",&mt,&mx)!=2) 
  { fprintf(stderr,"find-min: no entries\n");
    exit(1);
  }
  
  while (scanf("%lf%lf",&t,&x)==2)
  { if (x<mx) 
    { mx = x;
      mt = t;
    }
  }

  printf("%d\n",(int)(mt+0.5));

  exit(0);
}
