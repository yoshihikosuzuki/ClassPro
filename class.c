#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>

#include "libfastk.h"
#include "DB.h"
#include "bessel.h"
#include "ClassPro.h"

static const int    MAX_N_HC              = 5;
static const int    MAX_NITER             = 10;
static int          MIN_CNT_CHANGE        = 2;       // Every count change at a wall must be > this
static int          MAX_CNT_CHANGE        = 5;       // Every count change > this becomes a wall candidate   // TODO: change to sqrt(H-cov?)
static const double N_SIGMA_R             = 3;
static const double N_SIGMA_R_U           = 3;
static const int    N_BASE_EST            = 1000;
static const int    N_BASE_EST_MIN        = 1;
static const int    N_INTVL_EST           = 5;
static const int    N_INTVL_EST_MIN       = 1;
static const double MIN_P_NORMAL          = 0.1;     // min percentage of bases in reliable H/D intervals

static const double pethres_init[N_ETYPE] = {0.001, 0.05};
static const double pethres[N_ETYPE]      = {1e-10, 1e-5};

static const int    GAMMA_WEIGHT[2]       = {1, 1};   // Hyperparameters
static const int    ALPHA_WEIGHT          = 1;
static int          a_prior[2], b_prior[2], alpha_prior[2], lambda_prior[2];
static double       dg_sum_alpha_prior, eta_weight_k_prior[2], eta_const_k_prior[2];

static double       logfact[32768];

#define DIGAMMA_MAX 10000000   // TODO: more efficient way (doubling?)
static double       digamma[DIGAMMA_MAX+1];   // Only for integers

void precompute_probs()
{ for (int n = 1; n < 32768; n++)
    logfact[n] = logfact[n-1]+log(n);

  digamma[1] = -0.57721566490153;
  for (int n = 1; n < DIGAMMA_MAX; n++)
    digamma[n+1] = digamma[n]+(1./(double)n);

  return;
}

/*******************************************************************************************
 *
 *  SEQUENCE CONTEXT
 *
 ********************************************************************************************/

void calc_seq_context(Seq_Ctx *lctx, Seq_Ctx *rctx, char *seq, const int rlen)
{ int       in_hp, in_ds, in_ts;
  const int rlenm1 = rlen-1;

  in_ds = in_ts = 0;
  for (int i = 1; i < rlen; i++)
    { in_hp = (seq[i-1] == seq[i]) ? 1 : 0;
      in_ds = in_ts = 0;

      if (in_hp)
        { lctx[i][HP] = lctx[i-1][HP]+1;
          lctx[i][DS] = rctx[i-1][DS] = 0;
        }
      else
        { lctx[i][HP] = 1;
          lctx[i][DS] = rctx[i-1][DS] = 1;
          for (int j = i-lctx[i-1][HP], n = 0; j < i; j++, n++)
            rctx[j][HP] = lctx[i-1-n][HP];
          if (i >= 3 && seq[i-3] == seq[i-1] && seq[i-2] == seq[i])
            { lctx[i][DS] = lctx[i-2][DS]+1;
              in_ds = 1;
            }
        }

      if (!in_ds)
        { int l = i-1;
          while (lctx[l][DS] > 1)
            l--;
          if (l < i-1)
            for (int j = l-1, n = 0; j < i; j++, n++)
              rctx[j-1][DS] = lctx[i-1-n][DS];
        }

      if (i >= 2)
        { if (in_hp && seq[i-2] == seq[i-1])
            lctx[i][TS] = rctx[i-2][TS] = 0;
          else if (i >= 5 && seq[i-5] == seq[i-2] && seq[i-4] == seq[i-1] && seq[i-3] == seq[i])
            { lctx[i][TS] = lctx[i-3][TS]+1;
              in_ts = 1;
            }
          else
            lctx[i][TS] = rctx[i-1][TS] = rctx[i-2][TS] = 1;

          if (!in_ts)
            { int l = i-1;
              while (lctx[l][TS] > 1)
                l--;
              if (l < i-1)
                for (int j = l-2, n = 0; j < i; j++, n++)
                  rctx[j-2][TS] = lctx[i-1-n][TS];
            }
        }
    }

  for (int j = rlen-lctx[rlenm1][HP], n = 0; j < rlen; j++, n++)
    rctx[j][HP] = lctx[rlenm1-n][HP];

  if (in_ds)
    { int l = rlenm1;
      while (lctx[l][DS] > 1)
        l--;
      if (l < rlenm1)
        for (int j = l-1, n = 0; j < rlen; j++, n++)
          rctx[j-1][DS] = lctx[rlenm1-n][DS];
    }

  if (in_ts)
    { int l = rlenm1;
      while (lctx[l][TS] > 1)
        l--;
      if (l < rlenm1)
        for (int j = l-2, n = 0; j < rlen; j++, n++)
          rctx[j-2][TS] = lctx[rlenm1-n][TS];
    }

  rctx[rlenm1][DS] = rctx[rlenm1][TS] = rctx[rlen-2][TS] = 0;

#ifdef DEBUG_CTX
  const Seq_Ctx *ctx[N_WTYPE]  = {lctx, rctx};
  const char    *dir[N_WTYPE]  = {"L", "R"};
  const char    *name[N_CTYPE] = {"HP", "DS", "TS"};
  const int      W = 50;

  fprintf(stderr,"Seq  %.*s...%.*s\n",W,seq,W,seq+rlen-W);
  for (int t = HP; t <= TS; t++)
    { for (int d = DROP; d <= GAIN; d++)
        { fprintf(stderr,"%s %s ",name[t],dir[d]);
          for (int i = 0; i < W; i++)
            fprintf(stderr,"%u",ctx[d][i][t]);
          fprintf(stderr,"...");
          for (int i = rlen-W; i < rlen; i++)
            fprintf(stderr,"%u",ctx[d][i][t]);
          fprintf(stderr,"\n");
        }
    }
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  return;
}

/*******************************************************************************************
 *
 *  PROBABILITIES
 *
 ********************************************************************************************/

static inline double logp_poisson(int k, int lambda)
{ return k * log((double)lambda) - lambda - logfact[k];
}

static inline double logp_skellam(int k, double lambda)
{ return -2. * lambda + log(bessi(abs(k),2.*lambda));
}

static inline double logp_binom(int k, int n, double p)
{ return logfact[n] - logfact[k] - logfact[n-k] + k * log(p) + (n-k) * log(1-p);
}

static inline double logp_binom_pre(int k, int n, double lpe, double l1mpe)
{ return logfact[n] - logfact[k] - logfact[n-k] + k * lpe + (n-k) * l1mpe;
}

// TODO: chi-square when k,n are large
static inline double binom_test_g(int k, int n, double pe, int exact)
{ const double lpe   = log(pe);   // TODO: precompute in error models?
  const double l1mpe = log(1-pe);
  const double mean  = n*pe;
  const int decrease = ((double)k >= mean) ? 1 : 0;
  double p, p_first, p_curr;

/*#ifdef DEBUG_BINOM
  fprintf(stderr,"BinomTest(%d; %d, %lf) [%s]\n",k,n,pe,exact ? "exact" : "approx");
  fprintf(stderr,"    k (%d) %s mean (%lf)\n",k,decrease ? ">=" : "<",mean);
#endif*/

  if (decrease)
    { p = p_first = exp(logp_binom_pre(k,n,lpe,l1mpe));
      for (int x = k+1; x <= n; x++)
        { p += p_curr = exp(logp_binom_pre(x,n,lpe,l1mpe));
          if (exact == 0 && 10 * p_curr < p_first)
            break;
        }
    }
  else
    { p = p_first = (k == 0) ? 0. : exp(logp_binom_pre(k-1,n,lpe,l1mpe));
      for (int x = k-2; x >= 0; x--)
        { p += p_curr = exp(logp_binom_pre(x,n,lpe,l1mpe));
          if (exact == 0 && 10 * p_curr < p_first)
            break;
        }
      p = 1-p;
    }

  return p;
}

/*******************************************************************************************
 *
 *  ERROR MODEL PER POSITION
 *
 ********************************************************************************************/

static int CMAX;
#define CIDX(cout, cin) ((cout)*(((cout)+1))>>1)+(cin)-1

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

/*******************************************************************************************
 *
 *  WALL DETECTION
 *
 ********************************************************************************************/

static void check_drop(int i, uint16 *profile, int plen, Seq_Ctx *ctx[2], Error_Model *emodel, P_Error *perror, Error_Intvl *eintvl[N_ETYPE], const int eidx, const int K)
{ double max_p[N_ETYPE] = {-1., -1.};
  int    max_j[N_ETYPE] = {-1, -1};

  const double _pe = emodel[0].pe[1];
  double pe, ps, po;
  int m, n, j;
  
  const int ipk = i+K-1;

#ifdef DEBUG_ERROR
  fprintf(stderr,"j >= %d + n - m\n",ipk);
#endif

  pe = _pe;

  // High-complexity errors
  m = 0;
  for (n = 0; n <= MAX_N_HC; n++)
    { j = ipk+n-m;
      if (j >= plen)
        break;
      
      if (n > 1)
        pe *= _pe;

      // TODO: pe is used NOT here BUT BinomTest... => precompute binom for pe, pe^2, ..., pe^5 (ok for pe?)
      // TODO: Make sure both spe[i] and spe[j] are computed with the same error type
      ps = perror[i][SELF][DROP] * perror[j][SELF][GAIN] * pe;
      po = perror[i][OTHERS][DROP] * perror[j][OTHERS][GAIN] * pe;
      // TODO: check p_diff (i.e. Skellam)?

      if (max_p[0] < ps)
        { max_p[0] = ps;
          max_j[0] = j;
        }
      if (max_p[1] < po)
        { max_p[1] = po;
          max_j[1] = j;
        }

      // TODO: In the case of errors in others, find max(pe * p_diff) and return only pe? (to find a pair with a similar count change, as no count change gives probability 1)

#ifdef DEBUG_ERROR
      fprintf(stderr,"  @ j = %d (m = %d, n = %d): %d -> %d\n",j,m,n,profile[j-1],profile[j]);
      fprintf(stderr,"    Ps (= %lf * %lf * %lf) = %lf",perror[i][SELF][DROP],perror[j][SELF][GAIN],pe,ps);
      if (ps > pethres[SELF])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
      fprintf(stderr,"    Po (= %lf * %lf * %lf) = %lf",perror[i][OTHERS][DROP],perror[j][OTHERS][GAIN],pe,po);
      if (po > pethres[OTHERS])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
#endif
    }

  // Low-complexity errors   // TODO: consider only max prob error type?
  for (int t = 0; t < 3; t++)
    { m = ctx[0][i][t];
      if (m == 0)
        continue;

      int flag_of = 0;
      n = 0;
      while (1)
        { int idx = i+(t+1)*(n+1);
          if (idx >= plen)
            { flag_of = 1;
              break;
            }
          if (ctx[0][idx][t] != m+n+1)
            break;
          n++;
        }
      
      if (flag_of)
        { ps = perror[i][SELF][DROP] * perror[i][SELF][DROP];
          po = perror[i][OTHERS][DROP] * perror[i][OTHERS][DROP];

#ifdef DEBUG
          fprintf(stderr,"Drop ctx overflow!\n");
#endif
          continue;  // TODO: special treatment instaed of skip?
        }
      else
        { m *= t+1;
          n *= t+1;

          j = ipk+n-m;
          if (j >= plen)
            continue;

          if (profile[j-1]>profile[j])
            { perror[j][SELF][GAIN] = 0.;
              perror[j][OTHERS][GAIN] = 1.;
            }
          else
            { perror[j][SELF][GAIN] = binom_test_g(profile[j-1],profile[j],emodel[t].pe[MIN(ctx[1][j][t],emodel[t].lmax)],0);
              perror[j][OTHERS][GAIN] = binom_test_g(profile[j]-profile[j-1],profile[j],emodel[t].pe[MIN(ctx[1][j][t],emodel[t].lmax)],0);
            }

          ps = perror[i][SELF][DROP] * perror[j][SELF][GAIN];
          po = perror[i][OTHERS][DROP] * perror[j][OTHERS][GAIN];

#ifdef DEBUG_ERROR
      fprintf(stderr,"  @ j = %d (t = %d, m = %d, n = %d): %d -> %d\n",j,t,m,n,profile[j-1],profile[j]);
      fprintf(stderr,"    Ps (= %lf * %lf) = %lf",perror[i][SELF][DROP],perror[j][SELF][GAIN],ps);
      if (ps > pethres[SELF])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
      fprintf(stderr,"    Po (= %lf * %lf) = %lf",perror[i][OTHERS][DROP],perror[j][OTHERS][GAIN],po);
      if (po > pethres[OTHERS])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
#endif
        }

      if (max_p[0] < ps)
        { max_p[0] = ps;
          max_j[0] = j;
        }
      if (max_p[1] < po)
        { max_p[1] = po;
          max_j[1] = j;
        }
    }

  for (int e = SELF; e <= OTHERS; e++)
    { eintvl[e][eidx].i = i;
      eintvl[e][eidx].j = max_j[e];
      eintvl[e][eidx].pe = max_p[e];
    }

  return;
}

static void check_gain(int i, uint16 *profile, Seq_Ctx *ctx[2], Error_Model *emodel, P_Error *perror, Error_Intvl *eintvl[N_ETYPE], const int eidx, const int K)
{ double max_p[N_ETYPE] = {-1., -1.};
  int    max_j[N_ETYPE] = {-1, -1};

  const double _pe = emodel[0].pe[1];
  double pe, ps, po;
  int m, n, j;
  
  const int ipk = i-K+1;

#ifdef DEBUG_ERROR
  fprintf(stderr,"j <= %d - n + m\n",ipk);
#endif

  pe = _pe;

  // High-complexity errors
  m = 0;
  for (n = 0; n <= MAX_N_HC; n++)
    { j = ipk-n+m;
      if (j < 0)
        break;
      
      if (n > 1)
        pe *= _pe;

      // TODO: pe is used NOT here BUT BinomTest... => precompute binom for pe, pe^2, ..., pe^5 (ok for pe?)
      // TODO: Make sure both spe[i] and spe[j] are computed with the same error type
      ps = perror[j][SELF][DROP] * perror[i][SELF][GAIN] * pe;
      po = perror[j][OTHERS][DROP] * perror[i][OTHERS][GAIN] * pe;
      // TODO: check p_diff (i.e. Skellam)?
      // TODO: Multiply Pr(read start)^(# of count differences)?

      if (max_p[0] < ps)
        { max_p[0] = ps;
          max_j[0] = j;
        }
      if (max_p[1] < po)
        { max_p[1] = po;
          max_j[1] = j;
        }

#ifdef DEBUG_ERROR
      fprintf(stderr,"  @ j = %d (m = %d, n = %d): %d -> %d\n",j,m,n,profile[j-1],profile[j]);
      fprintf(stderr,"    Ps (= %lf * %lf * %lf) = %lf",perror[j][SELF][DROP],perror[i][SELF][GAIN],pe,ps);
      if (ps > pethres[SELF])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
      fprintf(stderr,"    Po (= %lf * %lf * %lf) = %lf",perror[j][OTHERS][DROP],perror[i][OTHERS][GAIN],pe,po);
      if (po > pethres[OTHERS])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
#endif
    }

  // Low-complexity errors   // TODO: consider only max prob error type?
  for (int t = 0; t < 3; t++)
    { m = ctx[1][i][t];
      if (m == 0)
        continue;

      int flag_of = 0;
      n = 0;
      while (1)
        { int idx = i-(t+1)*(n+1);
          if (idx < 0)
            { flag_of = 1;
              break;
            }
          if (ctx[1][idx][t] != m+n+1)
            break;
          n++;
        }
      
      if (flag_of)
        { ps = perror[i][SELF][GAIN] * perror[i][SELF][GAIN];
          po = perror[i][OTHERS][GAIN] * perror[i][OTHERS][GAIN];

#ifdef DEBUG
          fprintf(stderr,"Gain ctx overflow!\n");
#endif
          continue;   // TODO: special treatment instaed of skip?
        }
      else
        { m *= t+1;
          n *= t+1;

          j = ipk-n+m;
          if (j <= 0)
            continue;

          if (profile[j-1]<profile[j])
            { perror[j][SELF][DROP] = 0.;
              perror[j][OTHERS][DROP] = 1.;
            }
          else
            { perror[j][SELF][DROP] = binom_test_g(profile[j],profile[j-1],emodel[t].pe[MIN(ctx[0][j][t],emodel[t].lmax)],0);
              perror[j][OTHERS][DROP] = binom_test_g(profile[j-1]-profile[j],profile[j-1],emodel[t].pe[MIN(ctx[0][j][t],emodel[t].lmax)],0);
            }

          ps = perror[j][SELF][DROP] * perror[i][SELF][GAIN];
          po = perror[j][OTHERS][DROP] * perror[i][OTHERS][GAIN];

#ifdef DEBUG_ERROR
      fprintf(stderr,"  @ j = %d (t = %d, m = %d, n = %d): %d -> %d\n",j,t,m,n,profile[j-1],profile[j]);
      fprintf(stderr,"    Ps (= %lf * %lf) = %lf",perror[j][SELF][DROP],perror[i][SELF][GAIN],ps);
      if (ps > pethres[SELF])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
      fprintf(stderr,"    Po (= %lf * %lf) = %lf",perror[j][OTHERS][DROP],perror[i][OTHERS][GAIN],po);
      if (po > pethres[OTHERS])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
#endif
        }

      if (max_p[0] < ps)
        { max_p[0] = ps;
          max_j[0] = j;
        }
      if (max_p[1] < po)
        { max_p[1] = po;
          max_j[1] = j;
        }
    }

  for (int e = SELF; e <= OTHERS; e++)
    { eintvl[e][eidx].i = max_j[e];
      eintvl[e][eidx].j = i;
      eintvl[e][eidx].pe = max_p[e];
    }

  return;
}

static void correct_wall_cnt(int beg, int end, int *ci, int *cj, uint16 *profile, Seq_Ctx *ctx[N_WTYPE], const int K)
{ int n_gain = 0, n_drop = 0;
  int first, last, lmax;

  last = MIN(beg+K-1,end-1);
  for (int i = beg; i < last; i++)
    n_gain += MAX((int)profile[i+1]-profile[i],0);
  if (beg+K-1 < end)
    { lmax = 0;
      for (int t = HP; t <= TS; t++)
        { int l = ctx[GAIN][beg+K-1][t]*(t+1);
          if (lmax < l)
            lmax = l;
        }
      last = beg+lmax;
      for (int i = beg; i < last; i++)
        n_gain -= MAX((int)profile[i]-profile[i+1],0);
    }
  
  first = MAX(end-K+1,beg);
  for (int i = first; i < end-1; i++)
    n_drop += MAX((int)profile[i]-profile[i+1],0);
  if (beg < end-K+1)
    { lmax = 0;
      for (int t = HP; t <= TS; t++)
        { int l = ctx[DROP][end-K+1][t]*(t+1);
          if (lmax < l)
            lmax = l;
        }
      first = end-lmax;
      for (int i = first; i < end-1; i++)
        n_drop -= MAX((int)profile[i+1]-profile[i],0);
    }

  *ci = MIN(profile[beg]+MAX(n_gain,0),32767);
  *cj = MIN(profile[end-1]+MAX(n_drop,0),32767);

  last = MIN(beg+2*K,end);
  for (int i = beg; i < last; i++)
    if (*ci < profile[i])
      *ci = profile[i];
  first = MAX(end-2*K,beg);
  for (int i = first; i < end; i++)
    if (*cj < profile[i])
      *cj = profile[i];

#ifdef DEBUG_COR
  fprintf(stderr,"@ (%d,%d-1) %d(=%d+%d), %d(=%d+%d)\n",beg,end,*ci,profile[beg],n_gain,*cj,profile[end-1],n_drop);
#endif

  return;
}

void find_wall(uint16 *profile, int plen, Seq_Ctx *ctx[N_WTYPE], Error_Model *emodel, P_Error *perror, P_Error *cerror, Error_Intvl *eintvl[N_ETYPE], Intvl *intvl, Rel_Intvl *rintvl, int *wall, char *asgn, const int K, int *_N, int *_M)
{ int      N;
  int      cng, e, w;
  uint16   cin, cout;
  int      maxt, maxl;
  double   pe, maxpe;
  int      eidx, ridx;

  int THRES_R = 1000;   // TODO: calc from D-depth or seed selection

  for (int i = 0; i < plen; i++)
    { asgn[i] = 0;
      for (e = SELF; e <= OTHERS; e++)
        for (w = DROP; w <= GAIN; w++)
          perror[i][e][w] = cerror[i][e][w] = 0.;
    }

  eidx = 0;
  for (int i = 1; i < plen; i++)
    { w = (profile[i-1] >= profile[i]) ? DROP : GAIN;
      cng = abs((int)profile[i-1]-(int)profile[i]);

      if (w == DROP)
        { cin  = profile[i];
          cout = profile[i-1];
        }
      else
        { cin  = profile[i-1];
          cout = profile[i];
        }

      maxt = -1, maxl = -1;
      maxpe = 0.;
      for (int t = HP; t <= TS; t++)
        { int length = MIN(ctx[w][i][t],emodel[t].lmax);
          pe = emodel[t].pe[length];
          if (maxpe < pe)
            { maxpe = pe;
              maxt = t;
              maxl = length;
            }
        }

      int is_wall = 0;
      if (cng <= MIN_CNT_CHANGE)
        is_wall = 0;
      else if (cng > MAX_CNT_CHANGE)
        is_wall = 1;
      else if (cout <= CMAX && cin < MAX(emodel[maxt].cthres[maxl][cout][SELF],emodel[maxt].cthres[maxl][cout][OTHERS]))
        is_wall = 1;
      if (MIN(profile[i-1],profile[i]) >= THRES_R)
        is_wall = 0;

      if (is_wall == 0)   // Not a wall candidate
        continue;

      if (perror[i][SELF][w] == 0.)
        { //perror[i][SELF][w] = (cout <= CMAX) ? emodel[maxt].pe_bt[maxl][CIDX(cout,cin)] : binom_test_g(cin,cout,maxpe,0);
          perror[i][SELF][w] = binom_test_g(cin,cout,maxpe,0);
        }
      if (perror[i][OTHERS][w] == 0.)
        { //perror[i][OTHERS][w] = (cout <= CMAX) ? emodel[maxt].pe_bt[maxl][CIDX(cout,cout-cin)] : binom_test_g(cout-cin,cout,maxpe,0);
          perror[i][OTHERS][w] = binom_test_g(cout-cin,cout,maxpe,0);
        }

#ifdef DEBUG_ERROR
      const char _type = (cng == 0) ? '=' : ((cng > 0) ? '>' : '<');
      fprintf(stderr,"@ %d -> %d: %d -> %d (%c)",i-1,i,profile[i-1],profile[i],_type);
      
      for (int d = 0; d < 2; d++)
        { fprintf(stderr,"%s(",(d==0) ? "  L" : "  R");
          for (int t = 0; t < 3; t++)
            fprintf(stderr,"%d%c",ctx[d][i][t],(t<2) ? ',' : ')');
        }

      if (is_wall)
        { fprintf(stderr,"  (pe, ps, po) = (%lf, %lf, %lf)",maxpe,perror[i][SELF][w],perror[i][OTHERS][w]);
          if (asgn[i])
            fprintf(stderr," ***");
        }
      fprintf(stderr,"\n");
#endif

      if (perror[i][SELF][w] > pethres_init[SELF] || perror[i][OTHERS][w] < pethres_init[OTHERS])   // TODO: need this? init_thres are used in the count threshold precompute
        { asgn[i] = 1;

          if (w == DROP)
            { for (int n = 0; n <= MAX_N_HC; n++)
                { int j = i-1+K+n;
                  if (j >= plen || profile[j-1]>profile[j])
                    continue;
                  maxpe = 0.;
                  for (int t = 0; t < 3; t++)
                    { int length = MIN(ctx[1-w][j][t],emodel[t].lmax);
                      pe = emodel[t].pe[length];
                      if (maxpe < pe)
                        { maxpe = pe;
                          maxt = t;
                          maxl = length;
                        }
                    }
                  if (perror[j][SELF][1-w] == 0.)
                    { //perror[j][SELF][1-w] = (profile[j] <= CMAX) ? emodel[maxt].pe_bt[maxl][CIDX(profile[j],profile[j-1])] : binom_test_g(profile[j-1],profile[j],maxpe,0);
                      perror[j][SELF][1-w] = binom_test_g(profile[j-1],profile[j],maxpe,0);
                    }
                  if (perror[j][OTHERS][1-w] == 0.)
                    { //perror[j][OTHERS][1-w] = (profile[j] <= CMAX) ? emodel[maxt].pe_bt[maxl][CIDX(profile[j],profile[j]-profile[j-1])] : binom_test_g(profile[j]-profile[j-1],profile[j],maxpe,0);
                      perror[j][OTHERS][1-w] = binom_test_g(profile[j]-profile[j-1],profile[j],maxpe,0);
                    }
                }

              check_drop(i,profile,plen,ctx,emodel,perror,eintvl,eidx,K);
            }
          else
            { for (int n = 0; n <= MAX_N_HC; n++)
                { int j = i+1-K-n;
                  if (j <= 0 || profile[j-1]<profile[j])
                    continue;
                  maxpe = 0.;
                  for (int t = 0; t < 3; t++)
                    { int length = MIN(ctx[1-w][j][t],emodel[t].lmax);
                      pe = emodel[t].pe[length];
                      if (maxpe < pe)
                        { maxpe = pe;
                          maxt = t;
                          maxl = length;
                        }
                    }
                  if (perror[j][SELF][1-w] == 0.)
                    { //perror[j][SELF][1-w] = (profile[j-1] <= CMAX) ? emodel[maxt].pe_bt[maxl][CIDX(profile[j-1],profile[j])] : binom_test_g(profile[j],profile[j-1],maxpe,0);
                      perror[j][SELF][1-w] = binom_test_g(profile[j],profile[j-1],maxpe,0);
                    }
                  if (perror[j][OTHERS][1-w] == 0.)
                    { //perror[j][OTHERS][1-w] = (profile[j-1] <= CMAX) ? emodel[maxt].pe_bt[maxl][CIDX(profile[j-1],profile[j-1]-profile[j])] : binom_test_g(profile[j-1]-profile[j],profile[j-1],maxpe,0);
                      perror[j][OTHERS][1-w] = binom_test_g(profile[j-1]-profile[j],profile[j-1],maxpe,0);
                    }
                }

              check_gain(i,profile,ctx,emodel,perror,eintvl,eidx,K);
            }
          eidx++;
        }
    }

#ifdef DEBUG_INTVL
  const char etoc[2] = {'S', 'O'};
  fprintf(stderr,"# of intervals inspected = %d\n",eidx);
  for (int i = 0; i < eidx; i++)
    { for (int e = SELF; e <= OTHERS; e++)
        { double _pe = eintvl[e][i].pe;
          if (_pe < pethres[e])
            continue;
          
          int beg = eintvl[e][i].i;
          int end = eintvl[e][i].j;

          fprintf(stderr,"%c (%d,%d) (",etoc[e],beg,end);
          for (int t = HP; t <= TS; t++)
            fprintf(stderr,"%d%c",ctx[DROP][beg][t],(t<2) ? ',' : ')');
          fprintf(stderr,"-(");
          for (int t = HP; t <= TS; t++)
            fprintf(stderr,"%d%c",ctx[GAIN][end][t],(t<2) ? ',' : ')');
          fprintf(stderr," [%d|%d->%d|%d] %lf   \n",profile[beg-1],profile[beg],profile[end-1],profile[end],_pe);
        }
    }
#endif

  // TODO: use wall from the beginning?
  // Repeat counts cannot be walls
  THRES_R = lambda_prior[1] + 6*(int)sqrt(lambda_prior[1]);
  for (int i = 1; i < plen; i++)
    if (asgn[i] == 1 && MIN(profile[i-1],profile[i]) >= THRES_R)
      asgn[i] = 0;
  // Exclude walls that can be explained by errors in others
  e = OTHERS;
  for (int i = 0; i < eidx; i++)
    if (eintvl[e][i].pe >= pethres[e])
      asgn[eintvl[e][i].i] = asgn[eintvl[e][i].j] = 0;
  // Include walls that can be explained by errors in self
  e = SELF;
  for (int i = 0; i < eidx; i++)
    if (eintvl[e][i].pe >= pethres[e])
      asgn[eintvl[e][i].i] = asgn[eintvl[e][i].j] = 1;

  // Wall positions
  N = 0;   // # of intervals (N+1 is the length of `wall`)
  wall[N++] = 0;
  for (int i = 1; i < plen; i++)
    if (asgn[i] == 1)
      wall[N++] = i;
  wall[N] = plen;

#ifdef DEBUG_INTVL
  fprintf(stderr,"List of wall positions:\n    ");
  for (int i = 1; i < N; i++)
    fprintf(stderr,"%d ",wall[i]);
  fprintf(stderr,"\n");
#endif

  // Mask error positions
  for (int i = 1; i < N; i++)
    asgn[wall[i]] = 0;
  e = SELF;
  for (int i = 0; i < eidx; i++)
    if (eintvl[e][i].pe >= pethres[e])
      for (int j = eintvl[e][i].i; j < eintvl[e][i].j; j++)
        asgn[j] = 1;
  
  // Find reliable intervals and correct counts
  ridx = 0;
  int iidx = 0;
  for (int i = 0; i < N; i++)
    { int beg = wall[i];
      int end = wall[i+1];
      bool is_rel = false;
      bool is_err = false;

      // Reliable intervals must be long and not explained by errors in others
      if (asgn[beg] == 1 || asgn[end-1] == 1)
        { is_err = true;
        }
      else if (end-beg >= K)
        { int ci, cj;
          correct_wall_cnt(beg,end,&ci,&cj,profile,ctx,K);
          double _lambda = (double)(ci+cj)/2*(end-beg)/20000;
          double logp_sk = logp_skellam(ci-cj,_lambda);

#ifdef DEBUG_INTVL
          fprintf(stderr,"(%d,%d) [%d, %d] %lf",beg,end,ci,cj,exp(logp_sk));
#endif

          // Transition between wall counts can be explained by sampling fluctuation
          if (logp_sk >= -9.21)   // 0.0001
            { 
#ifdef DEBUG_INTVL
              fprintf(stderr,"   Reliable");
#endif
              
              is_rel = true;
              rintvl[ridx].i = beg;   // TODO: change to `b` and `e`?
              rintvl[ridx].j = end;
              rintvl[ridx].ci = ci;
              rintvl[ridx].cj = cj;   // TODO: change to `cjm1`?
              ridx++;

              if (beg > 0)
                { if (profile[beg-1] == ci)
                    { cerror[beg][SELF][DROP] = cerror[beg][SELF][GAIN] = 0.;
                      cerror[beg][OTHERS][DROP] = cerror[beg][OTHERS][GAIN] = 1.;
                    }
                  else
                    { w = (profile[beg-1] > ci) ? DROP : GAIN;
                      if (w == DROP)
                        { cin  = ci;
                          cout = profile[beg-1];
                        }
                      else
                        { cin  = profile[beg-1];
                          cout = ci;
                        }
                      maxt = -1, maxl = -1;
                      maxpe = 0.;
                      for (int t = HP; t <= TS; t++)
                        { int length = MIN(ctx[w][beg][t],emodel[t].lmax);
                          pe = emodel[t].pe[length];
                          if (maxpe < pe)
                            { maxpe = pe;
                              maxt = t;
                              maxl = length;
                            }
                        }
                      cerror[beg][SELF][w] = binom_test_g(cin,cout,maxpe,0);
                      cerror[beg][OTHERS][w] = binom_test_g(cout-cin,cout,maxpe,0);
                    }
                }

              if (end < plen)
                { if (cj == profile[end])
                    { cerror[end][SELF][DROP] = cerror[end][SELF][GAIN] = 0.;
                      cerror[end][OTHERS][DROP] = cerror[end][OTHERS][GAIN] = 1.;
                    }
                  else
                    { w = (cj > profile[end]) ? DROP : GAIN;
                      if (w == DROP)
                        { cin  = profile[end];
                          cout = cj;
                        }
                      else
                        { cin  = cj;
                          cout = profile[end];
                        }
                      maxt = -1, maxl = -1;
                      maxpe = 0.;
                      for (int t = HP; t <= TS; t++)
                        { int length = MIN(ctx[w][end][t],emodel[t].lmax);
                          pe = emodel[t].pe[length];
                          if (maxpe < pe)
                            { maxpe = pe;
                              maxt = t;
                              maxl = length;
                            }
                        }
                      cerror[end][SELF][w] = binom_test_g(cin,cout,maxpe,0);
                      cerror[end][OTHERS][w] = binom_test_g(cout-cin,cout,maxpe,0);
                    }
                }
            }
#ifdef DEBUG_INTVL
              fprintf(stderr,"\n");
#endif
        }

      intvl[iidx].i = beg;
      intvl[iidx].j = end;
      intvl[iidx].is_rel = is_rel;
      intvl[iidx].is_err = is_err;
      intvl[iidx].asgn = (is_err) ? ERROR : N_STATE;
      iidx++;
    }

  *_N = N;
  *_M = ridx;
  return;
}

/*******************************************************************************************
 *
 *  COUNT HISTOGRAM
 *
 ********************************************************************************************/

void process_global_hist(char *FK_ROOT, int COVERAGE)
{ Histogram *H;
  int64     *hist;

  H = Load_Histogram(FK_ROOT);
  if (H == NULL)
    { fprintf(stderr,"%s: Cannot open %s.hist\n",Prog_Name,FK_ROOT);
      exit(1);
    }
  Modify_Histogram(H,H->low,H->high,0);
  hist = H->hist;

  // Find H/D peaks (= mean depths)
  if (COVERAGE != -1)
    { lambda_prior[1] = COVERAGE;
      lambda_prior[0] = COVERAGE >> 1;
      if (VERBOSE)
        fprintf(stderr,"Specified (H,D) peaks = (%d,%d)\n",lambda_prior[0],lambda_prior[1]);
    }
  else
    { int        maxcnt;               // Tallest peak; Assuming this is H or D
      int64      maxpk;
      int        lmaxcnt, rmaxcnt;     // l for around maxcnt/2, r for around maxcnt*2
      int64      lmaxpk, rmaxpk;       // Assuming taller peak out of l and r is D or H
      int        is_lpeak, is_rpeak;   // Is l/r a peak (rather than monotonously decreasing)?
      double     m, s;

      maxcnt = maxpk = 0;
      for (int i = MAX(2,H->low); i < MIN(1000,H->high); i++)   // NOTE: Assuming coverage at most 1000
        { if (hist[i-1] < hist[i] && hist[i] > hist[i+1] && maxpk < hist[i])
            { maxcnt = i;
              maxpk = hist[i];
            }
        }
      if (maxcnt < 10)
        { fprintf(stderr,"Could not find any peaks @ >= 10 in the histogram. Revise data or use -c.");
          exit(1);
        }
      else if (VERBOSE)
        fprintf(stderr,"Found tallest peak at %d (%lld)\n",maxcnt,maxpk);

      m = (double)maxcnt/2;
      s = sqrt(m);
      lmaxcnt = lmaxpk = is_lpeak = 0;
      for (int i = (int)round(m-s); i <= (int)round(m+s); i++)
        { if (lmaxpk < hist[i])
            { lmaxcnt = i;
              lmaxpk = hist[i];
              is_lpeak = (hist[i-1] < hist[i] && hist[i] > hist[i+1]) ? 1 : 0;
            }
        }

      m = (double)maxcnt*2;
      s = sqrt(m);
      rmaxcnt = rmaxpk = is_rpeak = 0;
      for (int i = (int)round(m-s); i <= (int)round(m+s); i++)
        { if (rmaxpk < hist[i])
            { rmaxcnt = i;
              rmaxpk = hist[i];
              is_rpeak = (hist[i-1] < hist[i] && hist[i] > hist[i+1]) ? 1 : 0;
            }
        }

      if (lmaxpk > rmaxpk)   // maxcnt is D-peak and lmaxcnt is H-peak
        { lambda_prior[1] = maxcnt;
          lambda_prior[0] = is_lpeak ? lmaxcnt : (maxcnt >> 1);
        }
      else   // maxcnt is H-peak and rmaxcnt is D-peak
        { lambda_prior[0] = maxcnt;
          lambda_prior[1] = is_rpeak ? rmaxcnt : (maxcnt << 1);
        }

      if (VERBOSE)
        fprintf(stderr,"Estimated (H,D) peaks = (%d,%d)\n",lambda_prior[0],lambda_prior[1]);
    }

  // Determine hyperparameters
  double totpk[2];   // Estimated total k-mer frequency in H/D components
  double p;

  for (int k = 0; k < 2; k++)
    { b_prior[k] = GAMMA_WEIGHT[k];
      a_prior[k] = lambda_prior[k] * b_prior[k];

      p = exp(logp_poisson(lambda_prior[k],lambda_prior[k]));
      totpk[k] = (double)hist[lambda_prior[k]]/p;
    }
  const int minidx = (totpk[0] < totpk[1]) ? 0 : 1;
  alpha_prior[minidx] = ALPHA_WEIGHT;
  alpha_prior[1-minidx] = (int)(ALPHA_WEIGHT * (totpk[1-minidx] / totpk[minidx]));

  // Precompute values used in VI of PMM
  dg_sum_alpha_prior = digamma[alpha_prior[0] + alpha_prior[1]];
  for (int k = 0; k < 2; k++)
    { eta_weight_k_prior[k] = digamma[a_prior[k]] - log(b_prior[k]);
      eta_const_k_prior[k] = digamma[alpha_prior[k]] - dg_sum_alpha_prior - a_prior[k] / b_prior[k];
    }

  if (VERBOSE)
    { fprintf(stderr,"Hyperparameters etc. for PMM:\n");
      fprintf(stderr,"    lambda_prior       = (%d,%d)\n",lambda_prior[0],lambda_prior[1]);
      fprintf(stderr,"    a_prior            = (%d,%d)\n",a_prior[0],a_prior[1]);
      fprintf(stderr,"    b_prior            = (%d,%d)\n",b_prior[0],b_prior[1]);
      fprintf(stderr,"    alpha_prior        = (%d,%d)\n",alpha_prior[0],alpha_prior[1]);
      fprintf(stderr,"    dg_sum_alpha_prior = %lf\n",dg_sum_alpha_prior);
      fprintf(stderr,"    eta_weight_k_prior = (%lf,%lf)\n",eta_weight_k_prior[0],eta_weight_k_prior[1]);
      fprintf(stderr,"    eta_const_k_prior  = (%lf,%lf)\n",eta_const_k_prior[0],eta_const_k_prior[1]);
    }

  Free_Histogram(H);
}
 
int pmm_vi(uint16 *profile, uint16 *nprofile, int plen, double *eta, double lambda[2])
{ int    N;   // # of normal counts, = length of `nprofile` and `eta`/2
  int    is_converged;
  double a[2], b[2], alpha[2], eta_weight_k[2], eta_const_k[2];
  double eta_sum, dg_sum_alpha;

  const int ethres = lambda_prior[0]-3*(int)sqrt(lambda_prior[0]);
  const int rthres = lambda_prior[1]+3*(int)sqrt(lambda_prior[1]);

  for (int k = 0; k < 2; k++)
    { a[k]      = a_prior[k];
      b[k]      = b_prior[k];
      alpha[k]  = alpha_prior[k];
      lambda[k] = lambda_prior[k];
    }

  N = 0;
  for (int i = 0; i < plen; i++)
    if (ethres <= profile[i] && profile[i] <= rthres)
      nprofile[N++] = profile[i];
  if (N < 2)
    return N;

  for (int t = 0; t < MAX_NITER; t++)
    {
#ifdef DEBUG_PMM
      fprintf(stderr,"[PMM epoch %d]\n",t+1);
      fflush(stderr);
#endif

      // Update eta (assignment weights)
      if (t == 0)
        { dg_sum_alpha = dg_sum_alpha_prior;
          for (int k = 0; k < 2; k++)
            { eta_weight_k[k] = eta_weight_k_prior[k];
              eta_const_k[k] = eta_const_k_prior[k];
            }
        }
      else
        { if (alpha[0]+alpha[1] > DIGAMMA_MAX)
            { fprintf(stderr,"alpha[0]+alpha[1] (%lf) exceeded DIGAMMA_MAX\n",alpha[0]+alpha[1]);
              exit(1);
            }
          dg_sum_alpha = digamma[(int)(alpha[0] + alpha[1])];
          for (int k = 0; k < 2; k++)
            { if (a[k] > DIGAMMA_MAX)
                { fprintf(stderr,"a[%d] (%lf) exceeded DIGAMMA_MAX\n",k,a[k]);
                  exit(1);
                }
              eta_weight_k[k] = digamma[(int)a[k]] - log(b[k]);
              if (alpha[k] > DIGAMMA_MAX)
                { fprintf(stderr,"alpha[%d] (%lf) exceeded DIGAMMA_MAX\n",k,alpha[k]);
                  exit(1);
                }
              eta_const_k[k] = digamma[(int)alpha[k]] - dg_sum_alpha - a[k] / b[k];
            }
        }
      for (int n = 0; n < N; n++)
        { eta_sum = 0.;
          for (int k = 0; k < 2; k++)
            { int idx = (n<<1)|k;
              eta[idx] = exp(nprofile[n] * eta_weight_k[k] + eta_const_k[k]);
              eta_sum += eta[idx];
            }
          for (int k = 0; k < 2; k++)
            eta[(n<<1)|k] /= eta_sum;
        }

/*#ifdef DEBUG_PMM
      for (int n = 0; n < N; n++)
        fprintf(stderr,"n=%d, eta=(%lf,%lf)\n",n,eta[(n<<1)|0],eta[(n<<1)|1]);
#endif*/

      // Update distribution parameters
      // NOTE: Using expected value (as a good proxity for sum of assignments) in alpha updates
      for (int k = 0; k < 2; k++)
        a[k] = b[k] = 0.;
      for (int n = 0; n < N; n++)
        { for (int k = 0; k < 2; k++)
            { a[k] += eta[(n<<1)|k] * nprofile[n];
              b[k] += eta[(n<<1)|k];
            }
        }
      for (int k = 0; k < 2; k++)
        { alpha[k] = b[k] + alpha_prior[k];
          a[k] += a_prior[k];
          b[k] += b_prior[k];
        }

      // Compute lambda (depths)
      is_converged = 1;
      for (int k = 0; k < 2; k++)
        { double _lambda = a[k] / b[k];
          if (fabs(lambda[k] - _lambda) >= 0.1)
            is_converged = 0;
          lambda[k] = _lambda;
        }

#ifdef DEBUG_PMM
      fprintf(stderr,"  alpha  = (%lf,%lf)\n",alpha[0],alpha[1]);
      fprintf(stderr,"  a,b    = [(%lf,%lf), (%lf,%lf)]\n",a[0],b[0],a[1],b[1]);
      fprintf(stderr,"  lambda = (%lf,%lf)\n",lambda[0],lambda[1]);
      fflush(stderr);
#endif

      if (is_converged)
        break;
    }

#if defined(DEBUG_ITER) || defined(DEBUG_PMM)
  if (!is_converged)
    { fprintf(stderr,"PMM Not converged ");
#ifndef DEBUG_ITER
      fprintf(stderr,"\n");
#endif
      fflush(stderr);
    }
#endif

  // H-cov = D-cov / 2 if they are too close   // FIXME: H or D must be decided based on the global depth
  if (fabs(lambda[0]-lambda[1]) < sqrt(lambda[1]))
    { lambda[0] = lambda[1]/2;
#if defined(DEBUG_ITER) || defined(DEBUG_PMM)
      fprintf(stderr,"H=D/2 ");
#ifndef DEBUG_ITER
      fprintf(stderr,"\n");
#endif
      fflush(stderr);
#endif
    }

  return N;
}

/*******************************************************************************************
 *
 *  REFINE ASSIGNMENT
 *
 ********************************************************************************************/

static void nn_intvl(int idx, Rel_Intvl *rintvl, int M, char s, int ret[2])
{ int p = idx-1;
  while (p >= 0 && rintvl[p].asgn != s)
    p--;
  ret[0] = p;

  int n = idx+1;
  while (n < M && rintvl[n].asgn != s)
    n++;
  ret[1] = n;

  return;
}

static void nn_intvl_u(int idx, Intvl *intvl, int N, char s, int ret[2])
{ int p = idx-1;
  while (p >= 0 && intvl[p].asgn != s)
    p--;
  ret[0] = p;

  int n = idx+1;
  while (n < N && intvl[n].asgn != s)
    n++;
  ret[1] = n;

  return;
}

static void est_cnt_base(int idx, Intvl *intvl, int N, char s, uint16* profile, int ret[2])
{ int l, csum;

  l = N_BASE_EST;   // TODO: change to l=0 -> N_BASE_EST
  csum = 0;
  int p = idx-1;
  while (p >= 0 && l > 0) {
    if (intvl[p].asgn == s)
      for (int j = intvl[p].j-1; j >= intvl[p].i && l > 0; j--, l--)
        csum += profile[j];
    p--;
  }
  if (l > N_BASE_EST-N_BASE_EST_MIN)
    ret[0] = -1;
  else
    ret[0] = csum/(N_BASE_EST-l);
  
  l = N_BASE_EST;
  csum = 0;
  int n = idx+1;
  while (n < N && l > 0) {
    if (intvl[n].asgn == s)
      for (int j = intvl[n].i; j < intvl[n].j && l > 0; j++, l--)
        csum += profile[j];
    n++;
  }
  if (l > N_BASE_EST-N_BASE_EST_MIN)
    ret[1] = -1;
  else
    ret[1] = csum/(N_BASE_EST-l);

  return;
}

static void est_cnt_intvl(int idx, Rel_Intvl *rintvl, int M, char s, int ret[2])
{ int csum, lsum;
  int nadd, i;

  csum = lsum = 0;
  i = idx-1;
  nadd = 0;
  while (i >= 0 && nadd < N_INTVL_EST)
    { if (rintvl[i].asgn == s)
        { csum += (rintvl[i].j-rintvl[i].i)*(rintvl[i].ci+rintvl[i].cj)/2;
          lsum += rintvl[i].j-rintvl[i].i;
          nadd++;
        }
      i--;
    }
  if (nadd < N_INTVL_EST_MIN)
    ret[0] = -1;
  else
    ret[0] = csum/lsum;

  csum = lsum = 0;
  i = idx+1;
  nadd = 0;
  while (i < M && nadd < N_INTVL_EST)
    { if (rintvl[i].asgn == s)
        { csum += (rintvl[i].j-rintvl[i].i)*(rintvl[i].ci+rintvl[i].cj)/2;
          lsum += rintvl[i].j-rintvl[i].i;
          nadd++;
        }
      i++;
    }
  if (nadd < N_INTVL_EST_MIN)
    ret[1] = -1;
  else
    ret[1] = csum/lsum;

  return;
}

static void est_cnt_intvl_u(int idx, Intvl *intvl, int N, uint16 *profile, char s, int ret[2])
{ int csum, lsum;
  int nadd, i;

  csum = lsum = 0;
  i = idx-1;
  nadd = 0;
  while (i >= 0 && nadd < N_INTVL_EST)
    { if (intvl[i].asgn == s)
        { csum += (intvl[i].j-intvl[i].i)*(profile[intvl[i].i]+profile[intvl[i].j-1])/2;
          lsum += intvl[i].j-intvl[i].i;
          nadd++;
        }
      i--;
    }
  if (nadd < N_INTVL_EST_MIN)
    ret[0] = -1;
  else
    ret[0] = csum/lsum;

  csum = lsum = 0;
  i = idx+1;
  nadd = 0;
  while (i < N && nadd < N_INTVL_EST)
    { if (intvl[i].asgn == s)
        { csum += (intvl[i].j-intvl[i].i)*(profile[intvl[i].i]+profile[intvl[i].j-1])/2;
          lsum += intvl[i].j-intvl[i].i;
          nadd++;
        }
      i++;
    }
  if (nadd < N_INTVL_EST_MIN)
    ret[1] = -1;
  else
    ret[1] = csum/lsum;

  return;
}

static double calc_logp_e(int idx, Rel_Intvl *rintvl, int plen, P_Error *perror, int cov[])
{ Rel_Intvl ri = rintvl[idx];
  double logp_l, logp_r;
  double logp_po, logp_er;
  
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [ERROR]\n");
#endif

  logp_po = logp_poisson(ri.ci,cov[ERROR]);
  logp_er = -INFINITY;
  if (ri.i > 0)
    logp_er = log(perror[ri.i][SELF][DROP]);
  logp_l = MAX(logp_po,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [L] logp(PO)=%lf, logp(ER)=%lf\n",logp_po,logp_er);
#endif

  logp_po = logp_poisson(ri.cj,cov[ERROR]);
  logp_er = -INFINITY;
  if (ri.j < plen)
    logp_er = log(perror[ri.j][SELF][GAIN]);
  logp_r = MAX(logp_po,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [R] logp(PO)=%lf, logp(ER)=%lf\n",logp_po,logp_er);
#endif

  return logp_l+logp_r;
}

static double calc_logp_r(int idx, Rel_Intvl *rintvl, int M, int cov[])
{ Rel_Intvl ri = rintvl[idx];
  double logp_l, logp_r;
  double logp_sf, logp_er;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [REPEAT]\n");
#endif

  if (MAX(ri.ci,ri.cj) >= cov[REPEAT])
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,"    Larger than R-cov\n");
#endif
      return 0.;
    }

  int nn_idx[2];
  nn_intvl(idx,rintvl,M,DIPLO,nn_idx);
  int p = nn_idx[0];
  int n = nn_idx[1];
  int pe = -1, nb = -1;
  int pc = -1, nc = -1;

  if (p >= 0)
    { pe = rintvl[p].j;
      pc = rintvl[p].cj;
    }
  if (n < M)
    { nb = rintvl[n].i;
      nc = rintvl[n].ci;
    }
  if (p < 0 && n > M)
    { pc = nc = cov[DIPLO];
      pe = ri.i-10000;
      nc = ri.j+10000;
    }
  else if (p < 0)
    { pc = nc;
      pe = ri.i-(rintvl[n].i-ri.j);
    }
  else if (n > M)
    { nc = pc;
      nb = ri.j+(ri.i-rintvl[p].j);
    }

  double dr_ratio = 1.+N_SIGMA_R*(1./sqrt(cov[DIPLO]));
  pc = (int)((double)pc*dr_ratio);
  nc = (int)((double)nc*dr_ratio);

  double _lambda;
  _lambda = (double)pc*(ri.i-pe+1)/READ_LEN;
  logp_sf = logp_skellam(ri.ci-pc,_lambda);
  logp_er = (pc > ri.ci) ? logp_binom(ri.ci,pc,1-0.01) : -INFINITY;   // TODO: binom test? use ctx
  logp_l = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [L] logp(SF)=%lf, logp(ER)=%lf\n",logp_sf,logp_er);
#endif

  _lambda = (double)nc*(nb-ri.i+1)/READ_LEN;
  logp_sf = logp_skellam(nc-ri.cj,_lambda);
  logp_er = (ri.cj < nc) ? logp_binom(ri.cj,nc,1-0.01) : -INFINITY;
  logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [R] logp(SF)=%lf, logp(ER)=%lf\n",logp_sf,logp_er);
#endif

  return MAX(logp_l,logp_r);
}

static double calc_logp_hd(int s, int idx, Rel_Intvl *rintvl, int M, P_Error *cerror, int cov[])
{ Rel_Intvl ri = rintvl[idx];
  double logp_l = -INFINITY, logp_r = -INFINITY;
  double logp_sf, logp_er;

  int nn_idx[2];
  nn_intvl(idx,rintvl,M,s,nn_idx);
  int p = nn_idx[0];
  int n = nn_idx[1];

  if (p >= 0)
    { double _lambda = (double)cov[s]*(ri.i-rintvl[p].j+1)/READ_LEN;
      logp_sf = logp_skellam(ri.ci-rintvl[p].cj,_lambda);
      if (rintvl[p].j == ri.i)
        logp_er = log(MIN(cerror[ri.i][OTHERS][DROP],cerror[ri.i][OTHERS][GAIN]));   // TODO: MIN? MAX?
      else
        logp_er = -INFINITY;
      logp_l = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [L] logp(SF)=%lf, logp(ER)=%lf\n",logp_sf,logp_er);
#endif
    }
  if (n < M)
    { double _lambda = (double)cov[s]*(rintvl[n].i-ri.j+1)/READ_LEN;
      logp_sf = logp_skellam(rintvl[n].ci-ri.cj,_lambda);
      if (ri.j == rintvl[n].i)
        logp_er = log(MIN(cerror[ri.j][OTHERS][DROP],cerror[ri.j][OTHERS][GAIN]));
      else
        logp_er = -INFINITY;
      logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [R] logp(SF)=%lf, logp(ER)=%lf\n",logp_sf,logp_er);
#endif
    }

  if (p < 0 && n >= M)
    { logp_l = logp_poisson(ri.ci,cov[s]);
      logp_r = logp_poisson(ri.cj,cov[s]);
    }
  else if (p < 0)
    logp_l = logp_r;
  else if (n >= M)
    logp_r = logp_l;

  return MIN(logp_l,logp_r);
}

static double calc_logp_h(int idx, Rel_Intvl *rintvl, int M, P_Error *cerror, int cov[])
{ Rel_Intvl ri = rintvl[idx];

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [HAPLO]\n");
#endif

  int nn_idx[2];
  nn_intvl(idx,rintvl,M,DIPLO,nn_idx);
  int p = nn_idx[0];
  int n = nn_idx[1];
  if (p >= 0)
    { if (rintvl[p].cj < ri.ci)
        return -INFINITY;
    }
  if (n < M)
    { if (ri.cj > rintvl[n].ci)
        return -INFINITY;
    }
  
  int est_cnt[2];
  est_cnt_intvl(idx,rintvl,M,DIPLO,est_cnt);
  int pc = (int)((double)est_cnt[0]/1.25);   // TODO: change to N-sigma
  int nc = (int)((double)est_cnt[1]/1.25);

  if (pc > 0 && pc <= ri.ci)
    return -INFINITY;
  if (nc > 0 && nc <= ri.cj)
    return -INFINITY;

  return calc_logp_hd(HAPLO,idx,rintvl,M,cerror,cov);
}

static double calc_logp_d(int idx, Rel_Intvl *rintvl, int M, P_Error *cerror, int cov[])
{ Rel_Intvl ri = rintvl[idx];

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [DIPLO]\n");
#endif
 
  int est_cnt[2];
  est_cnt_intvl(idx,rintvl,M,HAPLO,est_cnt);
  int pc = (int)((double)est_cnt[0]*1.25);   // TODO: change to N-sigma
  int nc = (int)((double)est_cnt[1]*1.25);

  if (pc < 0 && nc < 0)
    pc = nc = (int)((double)cov[HAPLO]*1.25);
  else if (pc < 0)
   pc = nc;
  else if (nc < 0)
   nc = pc;
  
  if (ri.ci < pc && ri.cj < nc)
    return -INFINITY;

  return calc_logp_hd(DIPLO,idx,rintvl,M,cerror,cov);
}

static int update_state(int idx, Rel_Intvl *rintvl, int M, int plen, P_Error *perror, P_Error *cerror, int cov[])
{ char   s, smax = N_STATE;
  double logp, logpmax = -INFINITY;

  for (s = ERROR; s <= DIPLO; s++)
    { if (s == ERROR)
        logp = calc_logp_e(idx,rintvl,plen,perror,cov);
      else if (s == REPEAT)
        logp = calc_logp_r(idx,rintvl,M,cov);
      else if (s == HAPLO)
        logp = calc_logp_h(idx,rintvl,M,cerror,cov);
      else
        logp = calc_logp_d(idx,rintvl,M,cerror,cov);

      if (logp > logpmax)
        { smax = s;
          logpmax = logp;
        }

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,"idx=%d (%d,%d), s=%d, logp=%lf\n",idx,rintvl[idx].i,rintvl[idx].j,s,logp);
#endif
    }

#ifdef DEBUG
  if (smax > 4)
    { fprintf(stderr,"No valid probability for interval %d\n",idx);
      exit(1);
    }
#endif

  int changed = 0;
  if (rintvl[idx].asgn != smax)
    {
#ifdef DEBUG_REL
      fprintf(stderr,"state updated @ %d: %d -> %d\n",idx,rintvl[idx].asgn,smax);
#endif

      rintvl[idx].asgn = smax;
      changed = 1;
    }

  return changed;
}

typedef struct
  { int idx;
    int cnt;
  } Intvl_IC;

static int compare_iic(const void *a, const void *b)
{ return ((Intvl_IC*)a)->cnt - ((Intvl_IC*)b)->cnt;
}

void classify_reliable(Rel_Intvl *rintvl, int M, Intvl *intvl, int N, int plen, P_Error *perror, P_Error *cerror, int hcov, int dcov)
{ int cov[N_STATE] = {1, dcov+6*(int)sqrt(dcov), hcov, dcov};

  // Initial assignment
  for (int i = 0; i < M; i++)
    { int dmin = INT_MAX;
      for (int s = ERROR; s <= DIPLO; s++)
        { int d = abs(cov[s] - MAX(rintvl[i].ci,rintvl[i].cj));   // TODO: correct h/dcov estimates?
          if (dmin > d)
            { rintvl[i].asgn = s;
              dmin = d;
            }
        }
    }

#ifdef DEBUG_ITER
  char init_asgn[M];
  fprintf(stderr,"  Init : ");
  for (int i = 0; i < M; i++)
    { fprintf(stderr,"%c",stoc[(unsigned char)rintvl[i].asgn]);
      init_asgn[i] = rintvl[i].asgn;
    }
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  // Refine assignment
  Intvl_IC iord[M];
  for (int i = 0; i < M; i++)
    { iord[i].idx = i;
      iord[i].cnt = MIN(rintvl[i].ci,rintvl[i].cj);
    }
  qsort(iord,M,sizeof(Intvl_IC),compare_iic);

  int counter;
  for (counter = 0; counter < MAX_NITER; counter++)
    {
#ifdef DEBUG_REL
      fprintf(stderr,"[%dU]\n",counter+1);
#endif
      int changed = 0;
      for (int i = M - 1; i >= 0; i--)
        if (update_state(iord[i].idx,rintvl,M,plen,perror,cerror,cov))
          changed = 1;
      if (!changed) break;

#ifdef DEBUG_REL
      fprintf(stderr,"[%dD]\n",counter+1);
#endif

      changed = 0;
      for (int i = 0; i < M; i++)
        if (update_state(iord[i].idx,rintvl,M,plen,perror,cerror,cov))
          changed = 1;
      if (!changed) break;
    }

#ifdef DEBUG_ITER
  if (counter == MAX_NITER)
    fprintf(stderr,"REL Not converged\n");
  fprintf(stderr,"         ");
  for (int i = 0; i < M; i++)
    fprintf(stderr,"%c",(rintvl[i].asgn != init_asgn[i]) ? '*' : ' ');
  fprintf(stderr,"\n  Conv : ");
  for (int i = 0; i < M; i++)
    fprintf(stderr,"%c",stoc[(unsigned char)rintvl[i].asgn]);
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  // Copy assignments to `intvl`
  for (int ridx = 0, iidx = 0; ridx < M; ridx++, iidx++)
    { while (iidx < N && intvl[iidx].is_rel == 0)
        iidx++;
      if (iidx >= N || rintvl[ridx].i != intvl[iidx].i || rintvl[ridx].j != intvl[iidx].j)
        { fprintf(stderr,"Inconsistent reliable interval (%d,%d) != (%d,%d)\n",rintvl[ridx].i,rintvl[ridx].j,intvl[iidx].i,intvl[iidx].j);
          exit(1);
        }
      intvl[iidx].asgn = rintvl[ridx].asgn;
    }

  return;
}

static double calc_logp_e_u(int idx, Intvl *intvl, uint16 *profile, int plen, P_Error *perror, int cov[])
{ Intvl I = intvl[idx];
  double logp_l, logp_r;
  double logp_po, logp_er;
  
#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [ERROR]\n");
#endif

  logp_po = logp_poisson(profile[I.i],cov[ERROR]);
  logp_er = -INFINITY;
  if (I.i > 0)
    logp_er = log(perror[I.i][SELF][DROP]);
  logp_l = MAX(logp_po,logp_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [L] logp(PO)=%lf, logp(ER)=%lf\n",logp_po,logp_er);
#endif

  logp_po = logp_poisson(profile[I.j-1],cov[ERROR]);
  logp_er = -INFINITY;
  if (I.j < plen)
    logp_er = log(perror[I.j][SELF][GAIN]);
  logp_r = MAX(logp_po,logp_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [R] logp(PO)=%lf, logp(ER)=%lf\n",logp_po,logp_er);
#endif

  return logp_l+logp_r;
}

static double calc_logp_r_u(int idx, Intvl *intvl, int N, uint16 *profile, int cov[])
{ Intvl I = intvl[idx];
  double logp_l, logp_r;
  double logp_sf, logp_er;

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [REPEAT]\n");
#endif

  if (MAX(profile[I.i],profile[I.j-1]) >= cov[REPEAT])
    { 
#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
      fprintf(stderr,"    Larger than R-cov\n");
#endif
      return 0.;
    }

  int est_cnt[2];
  est_cnt_base(idx,intvl,N,DIPLO,profile,est_cnt);
  int pc = est_cnt[0];
  int nc = est_cnt[1];
  if (pc == -1 && nc == -1)
    { est_cnt_base(idx,intvl,N,HAPLO,profile,est_cnt);
      pc = est_cnt[0];
      nc = est_cnt[1];
      if (pc == -1 && nc == -1)
        pc = nc = cov[DIPLO];
      else if (pc == -1)
        pc = nc;
      else if (nc == -1)
        nc = pc;
    }
  else if (pc == -1)
    pc = nc;
  else if (nc == -1)
    nc = pc;

  double dr_ratio = 1.+N_SIGMA_R_U*(1./sqrt(cov[DIPLO]));
  pc = (int)((double)pc*dr_ratio);
  nc = (int)((double)nc*dr_ratio);

  if (pc <= profile[I.i] || profile[I.j-1] >= nc)
    return 0.;

  /*double _lambda;
  _lambda = (double)pc*(ri.i-pe+1)/READ_LEN;
  logp_sf = logp_skellam(ri.ci-pc,_lambda);*/
  logp_sf = -INFINITY;
  logp_er = (pc > profile[I.i]) ? logp_binom(profile[I.i],pc,1-0.01) : -INFINITY;   // TODO: binom test? use ctx
  logp_l = MAX(logp_sf,logp_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [L] logp(SF)=%lf, logp(ER)=%lf\n",logp_sf,logp_er);
#endif

  /*_lambda = (double)nc*(nb-ri.i+1)/READ_LEN;
  logp_sf = logp_skellam(nc-ri.cj,_lambda);*/
  logp_sf = -INFINITY;
  logp_er = (profile[I.j-1] < nc) ? logp_binom(profile[I.j-1],nc,1-0.01) : -INFINITY; 
  logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [R] logp(SF)=%lf, logp(ER)=%lf\n",logp_sf,logp_er);
#endif

  //return MAX(logp_l,logp_r);
  return logp_l+logp_r;
}

static double calc_logp_hd_u(int s, int idx, Intvl *intvl, int N, uint16 *profile, P_Error *perror, int cov[])
{ Intvl I = intvl[idx];
  double logp_l = -INFINITY, logp_r = -INFINITY;
  double logp_sf, logp_er;

  int nn_idx[2];
  nn_intvl_u(idx,intvl,N,s,nn_idx);
  int p = nn_idx[0];
  int n = nn_idx[1];

  if (p >= 0)
    { double _lambda = (double)cov[s]*(I.i-intvl[p].j+1)/READ_LEN;
      logp_sf = logp_skellam(profile[I.i]-profile[intvl[p].j-1],_lambda);
      if (intvl[p].j == I.i)
        logp_er = log(MIN(perror[I.i][OTHERS][DROP],perror[I.i][OTHERS][GAIN]));   // TODO: MIN? MAX?
      else
        logp_er = -INFINITY;
      logp_l = MAX(logp_sf,logp_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [L] logp(SF)=%lf, logp(ER)=%lf\n",logp_sf,logp_er);
#endif
    }
  if (n < N)
    { double _lambda = (double)cov[s]*(intvl[n].i-I.j+1)/READ_LEN;
      logp_sf = logp_skellam(profile[intvl[n].i]-profile[I.j-1],_lambda);
      if (I.j == intvl[n].i)
        logp_er = log(MIN(perror[I.j][OTHERS][DROP],perror[I.j][OTHERS][GAIN]));
      else
        logp_er = -INFINITY;
      logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"    [R] logp(SF)=%lf, logp(ER)=%lf\n",logp_sf,logp_er);
#endif
    }

  if (p < 0 || n >= N)
    { if (p < 0 && n >= N)   // TODO: allow short, isolated intervals to be H/D?
        { logp_l = logp_poisson(profile[I.i],cov[s]);
          logp_r = logp_poisson(profile[I.j-1],cov[s]);
        }
      else if (p < 0)
        logp_l = logp_r;
      else if (n >= N)
        logp_r = logp_l;
      return logp_l+logp_r;
    }

  //return MIN(logp_l,logp_r);
  return MAX(logp_l,logp_r);
}

static double calc_logp_h_u(int idx, Intvl *intvl, int N, uint16 *profile, P_Error *perror, int cov[])
{ Intvl I = intvl[idx];

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [HAPLO]\n");
#endif

  int est_cnt[2];
  est_cnt_intvl_u(idx,intvl,N,profile,DIPLO,est_cnt);
  int pc = (int)((double)est_cnt[0]/1.25);   // TODO: change to N-sigma
  int nc = (int)((double)est_cnt[1]/1.25);

  if (pc > 0 && pc <= profile[I.i])
    return -INFINITY;
  if (nc > 0 && nc <= profile[I.j-1])
    return -INFINITY;

  // TODO: NN H-interval + N-sigma check?

  return calc_logp_hd_u(HAPLO,idx,intvl,N,profile,perror,cov);
}

static double calc_logp_d_u(int idx, Intvl *intvl, int N, uint16 *profile, P_Error *perror, int cov[])
{ Intvl I = intvl[idx];

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [DIPLO]\n");
#endif
 
  int est_cnt[2];
  est_cnt_intvl_u(idx,intvl,N,profile,HAPLO,est_cnt);
  int pc = (int)((double)est_cnt[0]*1.25);   // TODO: change to N-sigma
  int nc = (int)((double)est_cnt[1]*1.25);

  if (pc < 0 && nc < 0)
    pc = nc = (int)((double)cov[HAPLO]*1.25);
  else if (pc < 0)
   pc = nc;
  else if (nc < 0)
   nc = pc;
  
  if (profile[I.i] < pc && profile[I.j-1] < nc)
    return -INFINITY;

  return calc_logp_hd_u(DIPLO,idx,intvl,N,profile,perror,cov);
}

static void update_state_u(int idx, Intvl *intvl, int N, uint16 *profile, int plen, P_Error *perror, int cov[])
{ char   s, smax = N_STATE;
  double logp, logpmax = -INFINITY;

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
      fprintf(stderr,"idx=%d (%d,%d)\n",idx,intvl[idx].i,intvl[idx].j);
#endif

  for (s = ERROR; s <= DIPLO; s++)
    { if (s == ERROR)
        logp = calc_logp_e_u(idx,intvl,profile,plen,perror,cov);
      else if (s == REPEAT)
        logp = calc_logp_r_u(idx,intvl,N,profile,cov);
      else if (s == HAPLO)
        logp = calc_logp_h_u(idx,intvl,N,profile,perror,cov);
      else
        logp = calc_logp_d_u(idx,intvl,N,profile,perror,cov);

      if (logp > logpmax)
        { smax = s;
          logpmax = logp;
        }

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
      fprintf(stderr,"idx=%d (%d,%d), s=%d, logp=%lf\n",idx,intvl[idx].i,intvl[idx].j,s,logp);
#endif
    }

#ifdef DEBUG
  if (smax >= N_STATE)
    { fprintf(stderr,"No valid probability for interval %d\n",idx);
      exit(1);
    }
#endif

  if (intvl[idx].asgn != smax)
    {
#ifdef DEBUG_UNREL
      fprintf(stderr,"state updated @ %d: %d -> %d\n",idx,intvl[idx].asgn,smax);
#endif

      intvl[idx].asgn = smax;
    }

  return;
}

void classify_unreliable(uint16 *profile, int plen, Intvl *intvl, int N, P_Error *perror, int hcov, int dcov)
{ int cov[N_STATE] = {1, dcov+6*(int)sqrt(dcov), hcov, dcov};

  // Check density of reliable H/D bases
  int rel_hd = 0;
  for (int i = 0; i < N; i++)
    { if (intvl[i].is_rel && (intvl[i].asgn == HAPLO || intvl[i].asgn == DIPLO))
        rel_hd += intvl[i].j - intvl[i].i;
    }
  double pnorm = (double)rel_hd/plen;
  if (pnorm < MIN_P_NORMAL)
    { 
#ifdef DEBUG_UNREL
      fprintf(stderr,"Too few reliable H/D bases.\n");
#endif
      
      for (int i = 0; i < N; i++)
        if (!intvl[i].is_rel && !intvl[i].is_err)
            intvl[i].asgn = REPEAT;

      return;
    }
  
  // Assignment
  Intvl_IC iord[N];
  for (int i = 0; i < N; i++)
    { iord[i].idx = i;
      iord[i].cnt = MIN(profile[intvl[i].i],profile[intvl[i].j-1]);
    }
  qsort(iord,N,sizeof(Intvl_IC),compare_iic);
  
  for (int i = N-1; i >= 0; i--)
    if (!intvl[iord[i].idx].is_rel && !intvl[iord[i].idx].is_err)
      update_state_u(iord[i].idx,intvl,N,profile,plen,perror,cov);

  for (int i = 0; i < N; i++)
    if (!intvl[iord[i].idx].is_rel && !intvl[iord[i].idx].is_err)
      update_state_u(iord[i].idx,intvl,N,profile,plen,perror,cov);

#ifdef DEBUG_ITER
  fprintf(stderr,"         ");
  for (int i = 0; i < N; i++)
    { if (!intvl[i].is_rel && !intvl[i].is_err)
        fprintf(stderr,"+");
      else
        fprintf(stderr," ");
    }
  fprintf(stderr,"\n  Intvl: ");
  for (int i = 0; i < N; i++)
    fprintf(stderr,"%c",stoc[(unsigned char)intvl[i].asgn]);
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  return;
}

// TODO: debug; make it work
void remove_slip(uint16 *profile, int plen, Seq_Ctx *ctx[N_WTYPE], char *crack)
{ int b, cnt, prev_d;
  bool in_hp;   // TODO: extend to DS/TS?
  char new_s;

  b = 0;
  cnt = 0;
  prev_d = -1;
  for (int i = 1; i < plen; i++)
    { int d = (profile[i-1] >= profile[i]) ? DROP : GAIN;
      if (d == DROP)
        in_hp = (ctx[d][i-1][HP] < ctx[d][i][HP]) ? true : false;
      else
        in_hp = (ctx[d][i-1][HP] > ctx[d][i][HP]) ? true : false;
      if (!in_hp || prev_d != d)
        { cnt = 0;
          continue;
        }
      if (crack[i-1] != crack[i])
        { cnt++;
          if (cnt == 1)
            b = i;
          else if (cnt >= 2)
            { if (d == DROP)
                new_s = (MIN(crack[i-1],crack[i]) == ERROR) ? ERROR : crack[i-1];
              else
                new_s = (MIN(crack[b-1],crack[b]) == ERROR) ? ERROR : crack[b];
#ifdef DEBUG_SLIP
              fprintf(stderr,"SLIP: [%d,%d) -> %c\n",b,i,stoc[(unsigned char)new_s]);
#endif
              for (int j = b; j < i; j++)
                crack[j] = new_s;
            }
        }
      prev_d = d;
    }

  return;
}

/*******************************************************************************************
 *
 *  SEED SELECTION
 *
 ********************************************************************************************/

/*static inline void find_seeds(uint16 *profile, char *asgn)
{
  return;
}*/
