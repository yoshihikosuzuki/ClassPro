#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "libfastk.h"
#include "DB.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static char *Usage = "[-s] [-e<int>] [-r<int>] <estimate>.class <truth>.class";

int main(int argc, char *argv[])
{ gzFile   fest, ftrue;
  kseq_t   *sest, *strue;

  bool SHOW_LQ = false, SHOW_CLASS = false;
  int THRES_LQ = -1;
  int THRES_R = 0;
  
  { int    i, j, k;
    int    flags[128];
    char  *eptr;
    (void) flags;

    ARG_INIT("class2acc");
    
    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
          { default:
              ARG_FLAGS("se")
              break;
            case 's':
              SHOW_CLASS = true;
              break;
            case 'e':
              SHOW_LQ = true;
              ARG_NON_NEGATIVE(THRES_LQ,"Min %%E-mer per read to show details")
              break;
            case 'r':
              ARG_NON_NEGATIVE(THRES_R,"Read with %%R-mer > this value is regarded as repeat")
              break;
          }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    if ((fest = gzopen(argv[1], "r")) == NULL)
      { fprintf(stderr,"%s: Cannot open %s [errno=%d]\n",Prog_Name,argv[1],errno);
        exit (1);
      }
    sest = kseq_init(fest);

    if ((ftrue = gzopen(argv[2], "r")) == NULL)
      { fprintf(stderr,"%s: Cannot open %s [errno=%d]\n",Prog_Name,argv[2],errno);
        exit (1);
      }
    strue = kseq_init(ftrue);
  }

  int id = 1;
  int64 ntot = 0, ncor = 0, nfne = 0;
  int64 ntot_normal = 0, ncor_normal = 0, nfne_normal = 0;
  int64 ntot_repeat = 0, ncor_repeat = 0, nfne_repeat = 0;
  int rcomp[4];
  while (kseq_read(sest) >= 0)
    { if (kseq_read(strue) < 0)
        { fprintf(stderr,"# seqs in %s > # seqs in %s\n",argv[1],argv[2]);
          exit(1);
        }

      if (strcmp(sest->name.s,strue->name.s) != 0)
        { fprintf(stderr,"Read %d inconsistent names: %s (estimate) vs %s (truth)\n",id,sest->name.s,strue->name.s);
          exit(1);
        }

      if (!(sest->seq.l == sest->qual.l && strue->seq.l == strue->qual.l && sest->seq.l == strue->seq.l))
        { fprintf(stderr,"Read %d inconsistent lengths\n",id);
          exit(1);
        }

      int i = 0;
      while (sest->qual.s[i] == 'N')
        { if (strue->qual.s[i] != 'N')
            { fprintf(stderr,"Read %d inconsistent # of prefix Ns (= K-1)\n",id);
              exit(1);
            }
          i++;
        }
      int rtot = strue->qual.l-i;
      int rcor = 0;
      int rfne = 0;
      for (int i = 0; i < 4; i++) rcomp[i] = 0;
      for (; i < (int)strue->qual.l; i++)
        { if (sest->qual.s[i] == strue->qual.s[i]) rcor++;
          if (strue->qual.s[i] == 'E' && sest->qual.s[i] != 'E') rfne++;
          switch (strue->qual.s[i])
            { case 'E':
                rcomp[0]++;
                break;
              case 'H':
                rcomp[1]++;
                break;
              case 'D':
                rcomp[2]++;
                break;
              case 'R':
                rcomp[3]++;
                break;
              default:
                fprintf(stderr,"Invalid class: %c\n",strue->qual.s[i]);
                break;
            }
        }
      ntot += rtot;
      ncor += rcor;
      nfne += rfne;

      if ((double)(rcomp[3])/rtot*100 > THRES_R)
        { ntot_repeat += rtot;
          ncor_repeat += rcor;
          nfne_repeat += rfne;
        }
      else
        { ntot_normal += rtot;
          ncor_normal += rcor;
          nfne_normal += rfne;
        }

      if (SHOW_LQ && (double)(rtot-rcor)/rtot*100 >= THRES_LQ)
        { fprintf(stdout,"Read %6d (%ld bp, %d classes): %%error = %4.1lf [%%E=%4.1lf,%%H=%4.1lf,%%D=%4.1lf,%%R=%4.1lf]\n",
                         id,strue->seq.l,rtot,(double)(rtot-rcor)/rtot*100,
                         (double)(rcomp[0])/rtot*100,(double)(rcomp[1])/rtot*100,(double)(rcomp[2])/rtot*100,(double)(rcomp[3])/rtot*100);
          if (SHOW_CLASS)
            fprintf(stdout,"  est: %s\ntruth: %s\n",sest->qual.s,strue->qual.s);
        }

      id++;
    }

  if (kseq_read(strue) >= 0)
    { fprintf(stderr,"# seqs in %s < # seqs in %s\n",argv[1],argv[2]);
      exit(1);
    }

  fprintf(stdout,"Accuracy = %4.2lf %% (= %lld / %lld), FN Error = %4.2lf %%\n",
                 (double)ncor/ntot*100,ncor,ntot,(double)nfne/ntot*100);
  fprintf(stdout,"[Normal] Accuracy = %4.2lf %% (= %lld / %lld), FN Error = %4.2lf %%\n",
                 (double)ncor_normal/ntot_normal*100,ncor_normal,ntot_normal,(double)nfne_normal/ntot_normal*100);
  fprintf(stdout,"[Repeat] Accuracy = %4.2lf %% (= %lld / %lld), FN Error = %4.2lf %%\n",
                 (double)ncor_repeat/ntot_repeat*100,ncor_repeat,ntot_repeat,(double)nfne_repeat/ntot_repeat*100);

  kseq_destroy(sest);
  gzclose(fest);
  kseq_destroy(strue);
  gzclose(ftrue);
  
  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  return 0;
}
