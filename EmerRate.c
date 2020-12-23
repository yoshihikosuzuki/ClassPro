#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#include "libfastk.h"

static char *Usage = "[-t<int(5)>] <source_root>[.prof]";

/****************************************************************************************
 *
 *  Test Stub
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ Profile_Index *P;
  int E_THRES;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("EmerRate");

    E_THRES = 5;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 't':
            ARG_POSITIVE(E_THRES,"E-mer count threshold")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc < 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  P = Open_Profiles(argv[1]);
  if (P == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  { int64   id;
    char   *eptr;
    uint16 *profile;
    int     plen, tlen;
    int     nemer;

    plen    = 20000;
    profile = Malloc(plen*sizeof(uint16),"Profile array");

    printf("RID\tRLEN\tN_EMER\tP_EMER\n");
    for (id = 1; id <= P->nreads; id++)
      { tlen = Fetch_Profile(P,id-1,plen,profile);
        if (tlen > plen)
          { plen    = 1.2*tlen + 1000;
            profile = Realloc(profile,plen*sizeof(uint16),"Profile array");
            Fetch_Profile(P,id-1,plen,profile);
          }
        nemer = 0;
        for (int i = 0; i < tlen; i++)
          if (profile[i] <= E_THRES) nemer++;
        printf("%lld\t%d\t%d\t%.3lf\n", id, tlen, nemer, (double)nemer/tlen);
      }
    free(profile);
  }

  Free_Profiles(P);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
