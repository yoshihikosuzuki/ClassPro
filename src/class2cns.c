#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "libfastk.h"
#include "DB.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static char *Usage = "<estimate>.class <fastk_root>[.prof]";

int main(int argc, char *argv[])
{ Profile_Index *P;
  int nthreads = 1;
  int nreads, nparts;
  int K, Km1;

  gzFile   fp;
  kseq_t   *seq;
  
  { int    i, j, k;
    int    flags[128];
    (void) flags;

    ARG_INIT("class2cns");
    
    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
          { default:
              ARG_FLAGS("")
              break;
          }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  { if ((fp = gzopen(argv[1], "r")) == NULL)
      { fprintf(stderr,"%s: Cannot open %s [errno=%d]\n",Prog_Name,argv[1],errno);
        exit (1);
      }
    seq = kseq_init(fp);
  }

  { if ((P = Open_Profiles(argv[2])) == NULL)
      { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,argv[2]);
        exit (1);
      }
    nreads = P->nreads;
    nparts = (nreads / nthreads) + (nreads % nthreads == 0 ? 0 : 1);
    K = P->kmer;
    Km1 = K-1;
  }

  while (kseq_read(seq) >= 0)
    { for (int i = Km1; i < (int)seq->seq.l; i++)
        { for (int j = i-Km1; j <= i; j++)
            printf("%c",seq->seq.s[j]);
          printf(" %c\n",seq->qual.s[i]);
        }
    }

  kseq_destroy(seq);
  gzclose(fp);
  
  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  return 0;
}
