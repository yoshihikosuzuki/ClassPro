#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "libfastk.h"
#include "DB.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define DEBUG

const char     stoc[4]    = { 'E', 'R', 'H', 'D' };
const char     ctos[128]        = { 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 3, 0, 0, 0,
                                    2, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 1, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0  };

static char *Usage = "\
[-s] [-e<int>] [-m<int(0)>] [-n<int(100)>] [-r<int(0)>] [-w<int>] [-p<read_profile[.prof]>] <estimate>.class <truth>.class\n\
\n\
  -e<int> : If specified with a value, classification information is shown for every read that has a misclassification rate larger than this value.\n\
\n\
  If `-e` is specified, then the following options become valid (Otherwise ignored):\n\
\n\
    -s           : If specified, for each read the ground-truth classification and the estimated classification are shown.\n\
    -m<int(0)>   : Minimum Repeat-mer rate of a read to be shown.\n\
    -n<int(100)> : Maximum Repeat-mer rate of a read to be shown.\n\
\n\
  -r<int(0)>   : Used for global accuracy calculation. Reads with a Repeat-mer rate larger than this value are regarded as 'Repeat reads'.\n\
  -w<int> : If specified with a value, instead of each read, accuracy is calculated for each window of the size of this value.\n\
  -p      : Path to .prof file.\n\
";

int main(int argc, char *argv[])
{ gzFile   fest, ftrue;
  kseq_t   *sest, *strue;
  Profile_Index *P = NULL;

  bool  SHOW_LQ    = false;   // -e
  bool  SHOW_CLASS = false;   // -s
  int   MIN_R      = 0;       // -m
  int   MAX_R      = 100;     // -n
  int   THRES_LQ   = -1;      // -e
  int   THRES_R    = 0;       // -r
  int   WINDOW     = -1;      // -w
  char *prof_root  = NULL;
  
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
            case 'm':
              ARG_NON_NEGATIVE(MIN_R,"Min %%R-mer per read to show details")
              break;
            case 'n':
              ARG_NON_NEGATIVE(MAX_R,"Max %%R-mer per read to show details")
              break;
            case 'e':
              SHOW_LQ = true;
              ARG_NON_NEGATIVE(THRES_LQ,"Min %%E-mer per read to show details")
              break;
            case 'r':
              ARG_NON_NEGATIVE(THRES_R,"Read with %%R-mer > this value is regarded as repeat")
              break;
            case 'w':
              ARG_NON_NEGATIVE(WINDOW,"Size of window = unit of coverage calculation")
              break;
            case 'p':
              prof_root = argv[i]+2;
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

    if (prof_root != NULL && (P = Open_Profiles(prof_root)) == NULL)
      { fprintf(stderr,"%s: Cannot open %s as a .prof file\n",Prog_Name,prof_root);
        exit (1);
      }
  }

  uint16 *profile = NULL;
  int     pmax = -1, plen = -1;
  int Km1 = -1;
  if (P != NULL)
    { Km1      = P->kmer-1;
      pmax     = 20000;
      profile  = Malloc(pmax*sizeof(uint16),"Profile array");
    }

  int id = 1;
  int64 ntot = 0, ncor = 0, nfne = 0;   // Total # of k-mers, correct classifications, false-negative error classifications
  int rtot, rcor, rfne;   // Per-read ...
  int wcor, wfne;         // Per-window ...
  int64 ntot_normal = 0, ncor_normal = 0, nfne_normal = 0;
  int64 ntot_repeat = 0, ncor_repeat = 0, nfne_repeat = 0;
  int rcomp[4];     // # of error/haplo/diplo/repeat-mers
  int wcomp[4];
  int64 cfm[4][4] = {};   // confusion matrix
  int64 scnts[2];   // Sum of counts of haplo/diplo-mers for coverage calculation
  double cov[2] = {-1, -1};    // per-read haplo/diplo coverages
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

      if (P != NULL)
        { plen = Fetch_Profile(P,(int64)id-1,pmax,profile);
          if (plen > pmax)
            { pmax    = 1.2*plen + 1000;
              profile = Realloc(profile,pmax*sizeof(uint16),"Profile array");
              Fetch_Profile(P,(int64)id-1,pmax,profile);
            }
          if (plen+Km1 != (int)sest->qual.l)
            { fprintf(stderr,"Read %d inconsist lengths: %ld (estimate) vs %d (profile)\n",id,sest->qual.l,plen+Km1);
              exit(1);
            }
        }

      int i = 0;
      while (sest->qual.s[i] == 'N')
        { if (strue->qual.s[i] != 'N')
            { fprintf(stderr,"Read %d inconsistent # of prefix Ns (= K-1)\n",id);
              exit(1);
            }
          i++;
        }
      rtot = strue->qual.l-i;
      rcor = rfne = 0;
      wcor = wfne = 0;
      for (int j = 0; j < 4; j++) 
        { rcomp[j] = 0;
          wcomp[j] = 0;
        }
      for (int j = 0; j < 2; j++) scnts[j] = 0;
      for (int c = 1; i < (int)strue->qual.l; i++, c++)
        { if (sest->qual.s[i] == strue->qual.s[i]) 
            { rcor++;
              wcor++;
            }
          if (strue->qual.s[i] == 'E' && sest->qual.s[i] != 'E')
            { rfne++;
              wfne++;
            }
          cfm[(int)ctos[(int)strue->qual.s[i]]][(int)ctos[(int)sest->qual.s[i]]]++;
          switch (strue->qual.s[i])
            { case 'E':
                rcomp[0]++; wcomp[0]++;
                break;
              case 'H':
                rcomp[1]++; wcomp[1]++;
                break;
              case 'D':
                rcomp[2]++; wcomp[2]++;
                break;
              case 'R':
                rcomp[3]++; wcomp[3]++;
                break;
              default:
                fprintf(stderr,"Invalid class: %c\n",strue->qual.s[i]);
                break;
            }
          if (P != NULL)
            { if (strue->qual.s[i] == 'H')
                scnts[0] += profile[i-Km1];
              else if (strue->qual.s[i] == 'D')
                scnts[1] += profile[i-Km1];
              if (WINDOW > 0)
                { if (c % WINDOW == 0)
                    { cov[0] = (wcomp[1] > 0) ? (double)(scnts[0])/wcomp[1] : -1;
                      cov[1] = (wcomp[2] > 0) ? (double)(scnts[1])/wcomp[2] : -1;
                      if (cov[0] == -1 || cov[1] == -1 || cov[0] > cov[1])
                        cov[0] = cov[1] = -1;
                      else 
                        cov[1] -= cov[0];
                      fprintf(stdout,"%%error = %4.1lf [H1-cov=%.lf,H2-cov=%.lf]\n",
                                     (double)(WINDOW-wcor)/WINDOW*100,cov[0],cov[1]);
                      for (int j = 0; j < 2; j++) scnts[j] = 0;
                      for (int j = 0; j < 4; j++) wcomp[j] = 0;
                      wcor = wfne = 0;
                    }
                }
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
      if (P != NULL)
        { cov[0] = (rcomp[1] > 0) ? (double)(scnts[0])/rcomp[1] : -1;
          cov[1] = (rcomp[2] > 0) ? (double)(scnts[1])/rcomp[2] : -1;
          if (cov[0] == -1 || cov[1] == -1 || cov[0] > cov[1])
            cov[0] = cov[1] = -1;
          else 
            cov[1] -= cov[0];
        }

      if (SHOW_LQ && (double)(rtot-rcor)/rtot*100 >= THRES_LQ
          && MIN_R <= (double)(rcomp[3])/rtot*100 && (double)(rcomp[3])/rtot*100 <= MAX_R)
        { fprintf(stdout,"Read %6d (%ld bp, %d classes): %%error = %4.1lf [%%E=%4.1lf,%%H=%4.1lf,%%D=%4.1lf,%%R=%4.1lf] [H1-cov=%.lf,H2-cov=%.lf]\n",
                         id,strue->seq.l,rtot,(double)(rtot-rcor)/rtot*100,
                         (double)(rcomp[0])/rtot*100,(double)(rcomp[1])/rtot*100,(double)(rcomp[2])/rtot*100,(double)(rcomp[3])/rtot*100,
                         cov[0],cov[1]);
          if (SHOW_CLASS)
            { fprintf(stdout,"truth: %s\n  est: ",strue->qual.s);
              for (int i = 0; i < (int)sest->qual.l; i++)
                if (strue->qual.s[i] != sest->qual.s[i])
                  fprintf(stdout,"%c",sest->qual.s[i]);
                else
                  fprintf(stdout,"-");
              fprintf(stdout,"\n");
            }
        }

      id++;
    }

  if (kseq_read(strue) >= 0)
    { fprintf(stderr,"# seqs in %s < # seqs in %s\n",argv[1],argv[2]);
      exit(1);
    }

  fprintf(stdout,"\nConfusion Matrix (Truth\\Est):\n  ");
  for (int i = 0; i < 4; i++)
    fprintf(stdout,"%15c",stoc[i]);
  fprintf(stdout,"\n");
  for (int i = 0; i < 4; i++)
    { fprintf(stdout,"%c:",stoc[i]);
      for (int j = 0; j < 4; j++)
        fprintf(stdout,"%15lld",cfm[i][j]);
      fprintf(stdout,"\n");
    }

  fprintf(stdout,"\nAccuracy = %4.2lf %% (= %lld / %lld), FN Error = %4.2lf %%\n",
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
