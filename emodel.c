#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "ClassPro.h"

int CMAX;

Error_Model *calc_init_thres()
{ Error_Model *emodel;
  int          cout, cin;
  double       pe, lpe, l1mpe, psum;
  bool         is_found[N_THRES][N_ETYPE];
  uint8        ct[N_ETYPE];

  CMAX = GLOBAL_COV[REPEAT];
  if (CMAX > 255)
    { fprintf(stderr,"Too high coverage\n");
      exit(1);
    }

  emodel = Malloc(sizeof(Error_Model)*N_CTYPE,"Allocating error model");

  for (int t = HP; t <= TS; t++)
    { emodel[t].lmax = (uint8)(MAX_N_LC/(t+1));

      // TODO: change to loading table
      emodel[t].pe = Malloc(sizeof(double)*(emodel[t].lmax+1),"Allocating pe");
      emodel[t].pe[0] = 0.;
      for (int l = 1; l <= emodel[t].lmax; l++)
        emodel[t].pe[l] = 0.002 * l * l + 0.002;

      emodel[t].cthres = Malloc(sizeof(uint8***)*(emodel[t].lmax+1),"Allocating cthres");
      for (int l = 1; l <= emodel[t].lmax; l++)
        { emodel[t].cthres[l] = Malloc(sizeof(uint8**)*CMAX,"Allocating cthres l");
          for (int c = 1; c < CMAX; c++)
            { emodel[t].cthres[l][c] = Malloc(sizeof(uint8*)*N_THRES,"Allocating cthres pe");
              for (int s = INIT; s <= FINAL; s++)
                emodel[t].cthres[l][c][s] = Malloc(sizeof(uint8)*N_ETYPE,"Allocating cthres thres");
            }
        }
    }

#ifdef DEBUG_EMODEL
  fprintf(stderr,"Thresholds for initial wall filtering (cmax = %d):\n",CMAX);
  fprintf(stderr,"          cout       :");
  for (cout = 1; cout < CMAX; cout++)
    fprintf(stderr," %3d",cout);
  fprintf(stderr,"\n  ( t, l, pe)\n");
#endif

  for (int t = HP; t <= TS; t++)
    for (int l = 1; l <= emodel[t].lmax; l++)
      { pe = emodel[t].pe[l];
        lpe = log(pe);
        l1mpe = log(1-pe);

#ifdef DEBUG_EMODEL
      fprintf(stderr,"  (%2d,%2d, %.3lf)\n",t,l,pe);
#endif

        for (cout = 1; cout < CMAX; cout++)
          { // initialize
            ct[SELF] = cout;
            ct[OTHERS] = 0;
            for (int s = INIT; s <= FINAL; s++)
              for (int e = SELF; e <= OTHERS; e++)
                { emodel[t].cthres[l][cout][s][e] = ct[e];
                  is_found[s][e] = false;
                }
            // find thresholds
            psum = 1.;
            for (cin = 0; cin <= cout; cin++)
              { if (is_found[INIT][SELF] && is_found[FINAL][SELF]
                    && is_found[INIT][OTHERS] && is_found[FINAL][OTHERS])
                  break;
                ct[SELF] = cin;
                ct[OTHERS] = cout-cin;
                psum -= exp(logp_binom_pre(cin,cout,lpe,l1mpe));
                for (int s = INIT; s <= FINAL; s++)
                  for (int e = SELF; e <= OTHERS; e++)
                    if (!is_found[s][e] && psum < PE_THRES[s][e])
                      { emodel[t].cthres[l][cout][s][e] = ct[e];
                        is_found[s][e] = true;
                      }
              }
          }

#ifdef DEBUG_EMODEL
        fprintf(stderr,"          cin(S_init ):");
        for (cout = 1; cout < CMAX; cout++)
          fprintf(stderr," %3d",emodel[t].cthres[l][cout][INIT][SELF]);
        fprintf(stderr,"\n          cin(S_final):");
        for (cout = 1; cout < CMAX; cout++)
          fprintf(stderr," %3d",emodel[t].cthres[l][cout][FINAL][SELF]);
        fprintf(stderr,"\n          cin(O_init ):");
        for (cout = 1; cout < CMAX; cout++)
          fprintf(stderr," %3d",emodel[t].cthres[l][cout][INIT][OTHERS]);
        fprintf(stderr,"\n          cin(O_final):");
        for (cout = 1; cout < CMAX; cout++)
          fprintf(stderr," %3d",emodel[t].cthres[l][cout][FINAL][OTHERS]);
        fprintf(stderr,"\n");
        fflush(stderr);
#endif
      }

  return emodel;
}

void free_emodel(Error_Model *emodel)
{ for (int t = HP; t <= TS; t++)
    { free(emodel[t].pe);
      for (int l = 1; l <= emodel[t].lmax; l++)
        { for (int c = 1; c < CMAX; c++)
            { for (int s = INIT; s <= FINAL; s++)
                free(emodel[t].cthres[l][c][s]);
              free(emodel[t].cthres[l][c]);
            }
          free(emodel[t].cthres[l]);
        }
      free(emodel[t].cthres);
    }
  free(emodel);

  return;
}
