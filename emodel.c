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
  double       pe, lpe, l1mpe;
  double       psum;

  const int MAX_CLEN = 20;
  CMAX = MIN(255,lambda_prior[1] + (int)(3*sqrt(lambda_prior[1])));
  const int cthres_size = (CMAX*(CMAX+3))>>1;

  // TODO: change to loading table
  emodel = Malloc(sizeof(Error_Model)*N_CTYPE,"Allocating error model");

  for (int t = HP; t <= TS; t++)
    { emodel[t].lmax = (uint8)(MAX_CLEN/(t+1));
      
      emodel[t].pe = Malloc(sizeof(double)*(emodel[t].lmax+1),"Allocating pe");
      emodel[t].pe[0] = 0.;
      for (int l = 1; l <= emodel[t].lmax; l++)
        emodel[t].pe[l] = 0.002 * l * l + 0.002;

      emodel[t].pe_bt = (double**)Malloc(sizeof(double*)*(emodel[t].lmax+1),"Allocating cthres");
      for (int l = 1; l <= emodel[t].lmax; l++)
        emodel[t].pe_bt[l] = (double*)Malloc(sizeof(double)*(cthres_size),"Allocating cthres pe");
      
      emodel[t].cthres = Malloc(sizeof(int**)*(emodel[t].lmax+1),"Allocating cthres");
      for (int l = 1; l <= emodel[t].lmax; l++)
        { emodel[t].cthres[l] = Malloc(sizeof(int*)*(CMAX+1),"Allocating cthres pe");
          for (int c = 1; c <= CMAX; c++)
            emodel[t].cthres[l][c] = Malloc(sizeof(int)*2,"Allocating cthres pe");
        }
    }

#ifdef DEBUG_BINOM
  fprintf(stderr,"Thresholds for initial wall filtering (cmax = %d; Given (c_out, c_in), it is a wall candidate if c_in <= this value) (cthres_size=%d):\n",CMAX,cthres_size);
  fprintf(stderr,"          cout  :");
  for (cout = 1; cout <= CMAX; cout++)
    fprintf(stderr," %3d",cout);
  fprintf(stderr,"\n  ( t, l, pe)\n");
#endif

  for (int t = HP; t <= TS; t++)
    for (int l = 1; l <= emodel[t].lmax; l++)
      { pe = emodel[t].pe[l];
        lpe = log(pe);
        l1mpe = log(1-pe);

#ifdef DEBUG_BINOM
      fprintf(stderr,"  (%2d,%2d, %lf)\n",t,l,pe);
#endif

        int idx = 0;
        for (cout = 1; cout <= CMAX; cout++)
          { emodel[t].cthres[l][cout][SELF] = cout;
            emodel[t].cthres[l][cout][OTHERS] = 0;

            int s_found = 0, o_found = 0;
            psum = 1.;
            for (cin = 0; cin <= cout; cin++)
              { emodel[t].pe_bt[l][idx] = (psum > 0.) ? psum : 0.;
                idx++;
                psum -= exp(logp_binom_pre(cin,cout,lpe,l1mpe));

                if (s_found == 0 && psum < pethres_init[SELF])
                  { emodel[t].cthres[l][cout][SELF] = cin;
                    s_found = 1;
                  }
                if (o_found == 0 && psum < pethres_init[OTHERS])
                  { emodel[t].cthres[l][cout][OTHERS] = cout - cin;
                    o_found = 1;
                  }
              }
          }

#ifdef DEBUG_BINOM
        fprintf(stderr,"\nlast idx=%d\n",idx);

        fprintf(stderr,"          cin(S):");
        for (cout = 1; cout <= CMAX; cout++)
          fprintf(stderr," %3d",emodel[t].cthres[l][cout][SELF]);
        fprintf(stderr,"\n          cin(O):");
        for (cout = 1; cout <= CMAX; cout++)
          fprintf(stderr," %3d",emodel[t].cthres[l][cout][OTHERS]);
        fprintf(stderr,"\n");
        for (cout = 1; cout <= CMAX; cout++)
          { for (cin = 0; cin <= cout; cin++)
              { int _idx = CIDX(cout,cin);
                double pe_ex = binom_test_g(cin,cout,emodel[t].pe[MIN(l,emodel[t].lmax)],0);
                fprintf(stderr,"(%d,%d)@%d=%lf(%lf) ",cout,cin,_idx,emodel[t].pe_bt[l][_idx],pe_ex);
              }
            fprintf(stderr,"\n");
          }
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
        { free(emodel[t].pe_bt[l]);
          for (int c = 1; c <= CMAX; c++)
            free(emodel[t].cthres[l][c]);
          free(emodel[t].cthres[l]);
        }
      free(emodel[t].pe_bt);
      free(emodel[t].cthres);
    }
  free(emodel);

  return;
}
