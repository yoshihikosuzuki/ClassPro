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
#include "ClassPro.h"

static const int    REL_MAX_NITER         = 10;
static const double N_SIGMA_R             = 3;
static const double N_SIGMA_R_U           = 3;
static const int    N_BASE_EST            = 1000;
static const int    N_BASE_EST_MIN        = 1;
static const int    N_INTVL_EST           = 5;
static const int    N_INTVL_EST_MIN       = 1;
static const double MIN_P_NORMAL          = 0.1;     // min percentage of bases in reliable H/D intervals

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
  fprintf(stderr,"  [E]");
#endif

  logp_po = logp_poisson(ri.ci,cov[ERROR]);
  logp_er = -INFINITY;
  if (ri.i > 0)
    logp_er = log(perror[ri.i][SELF][DROP]);
  logp_l = MAX(logp_po,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr," [L] logp(PO) = %.1lf, logp(ER) = %.1lf\n",logp_po,logp_er);
#endif

  logp_po = logp_poisson(ri.cj,cov[ERROR]);
  logp_er = -INFINITY;
  if (ri.j < plen)
    logp_er = log(perror[ri.j][SELF][GAIN]);
  logp_r = MAX(logp_po,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"      [R] logp(PO) = %.1lf, logp(ER) = %.1lf\n",logp_po,logp_er);
#endif

  return logp_l+logp_r;
}

static double calc_logp_r(int idx, Rel_Intvl *rintvl, int M, int cov[])
{ Rel_Intvl ri = rintvl[idx];
  double logp_l, logp_r;
  double logp_sf, logp_er;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [R]");
#endif

  if (MAX(ri.ci,ri.cj) >= cov[REPEAT])
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr," (Larger than R-cov)\n          logp(R) = 0.\n");
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
  fprintf(stderr," [L] logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
#endif

  _lambda = (double)nc*(nb-ri.i+1)/READ_LEN;
  logp_sf = logp_skellam(nc-ri.cj,_lambda);
  logp_er = (ri.cj < nc) ? logp_binom(ri.cj,nc,1-0.01) : -INFINITY;
  logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"      [R] logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
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
        logp_er = log(MAX(cerror[ri.i][OTHERS][DROP],cerror[ri.i][OTHERS][GAIN]));   // TODO: check this prob is what we expect
      else
        logp_er = -INFINITY;
      logp_l = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr," [L] logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
#endif
    }
  if (n < M)
    { double _lambda = (double)cov[s]*(rintvl[n].i-ri.j+1)/READ_LEN;
      logp_sf = logp_skellam(rintvl[n].ci-ri.cj,_lambda);
      if (ri.j == rintvl[n].i)
        logp_er = log(MAX(cerror[ri.j][OTHERS][DROP],cerror[ri.j][OTHERS][GAIN]));
      else
        logp_er = -INFINITY;
      logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"      [R] logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
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
  fprintf(stderr,"  [H]");
#endif

  // FIXME: This hard assignment must be inappropriate
  /*int nn_idx[2];
  nn_intvl(idx,rintvl,M,DIPLO,nn_idx);
  int p = nn_idx[0];
  int n = nn_idx[1];
  if (p >= 0)
    { if (rintvl[p].cj < ri.ci)
        { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
          fprintf(stderr," (Start > than left nearest D)\n          logp(H) = -inf\n");
#endif
          return -INFINITY;
        }
    }
  if (n < M)
    { if (ri.cj > rintvl[n].ci)
        { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
          fprintf(stderr," (End > than right nearest D)\n          logp(H) = -inf\n");
#endif
          return -INFINITY;
        }
    }
  
  // FIXME: This hard assignment must be inappropriate
  int est_cnt[2];
  est_cnt_intvl(idx,rintvl,M,DIPLO,est_cnt);
  int pc = (int)((double)est_cnt[0]/1.25);   // TODO: change to N-sigma
  int nc = (int)((double)est_cnt[1]/1.25);

  if (pc > 0 && pc <= ri.ci)
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr," (Start count > left est D / 1.25)\n          logp(H) = -inf\n");
#endif
      return -INFINITY;
    }
  if (nc > 0 && nc <= ri.cj)
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr," (End count > right est D / 1.25)\n          logp(H) = -inf\n");
#endif
      return -INFINITY;
    }*/

  return calc_logp_hd(HAPLO,idx,rintvl,M,cerror,cov);
}

static double calc_logp_d(int idx, Rel_Intvl *rintvl, int M, P_Error *cerror, int cov[])
{ Rel_Intvl ri = rintvl[idx];

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [D]");
#endif

  // FIXME: This hard assignment must be inappropriate
  /*int est_cnt[2];
  est_cnt_intvl(idx,rintvl,M,HAPLO,est_cnt);
  int pc = (int)((double)est_cnt[0]*1.25);   // TODO: change to N-sigma
  int nc = (int)((double)est_cnt[1]*1.25);

  if (pc < 0 && nc < 0)
    pc = nc = (int)((double)cov[HAPLO]*1.25);
  else if (pc < 0)
   pc = nc;
  else if (nc < 0)
   nc = pc;
  
  // FIXME: This hard assignment must be inappropriate
  if (ri.ci < pc && ri.cj < nc)
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr," (Start count < left est H*1.25 && End count < right est H*1.25)\n          logp(D) = -inf\n");
#endif
      return -INFINITY;
    }*/

  return calc_logp_hd(DIPLO,idx,rintvl,M,cerror,cov);
}

// TODO: array of functions
/*typedef double (*Func_Logp)(int idx, Rel_Intvl *rintvl, int plen, int M, P_Error *perror, P_Error *cerror, int cov[]);
Func_Logp calc_logps[N_STATE] = { &calc_logp_e, &calc_logp_r, &calc_logp_h, $calc_logp_d };*/

static inline bool update_state(int idx, Rel_Intvl *rintvl, int M, int plen, P_Error *perror, P_Error *cerror, int cov[], int counter)
{ int    s, smax = N_STATE;
  double logpmax = -INFINITY;
  double logps[N_STATE];

#ifdef DEBUG_REL
  { Rel_Intvl r    = rintvl[idx];
    const int iter = abs(counter);
    const char dir = (counter>0) ? 'U' : 'D';
    fprintf(stderr,"\n%d%c RI[%d] ",iter,dir,idx);
    fprintf(stderr,"@ %d -> %d: %d -> %d (%c)\n",r.i,r.j,r.ci,r.cj,stoc[(unsigned char)r.asgn]);
  }
#endif

  for (s = ERROR; s <= DIPLO; s++)
    { if (s == ERROR)
        logps[s] = calc_logp_e(idx,rintvl,plen,perror,cov);
      else if (s == REPEAT)
        logps[s] = calc_logp_r(idx,rintvl,M,cov);
      else if (s == HAPLO)
        logps[s] = calc_logp_h(idx,rintvl,M,cerror,cov);
      else
        logps[s] = calc_logp_d(idx,rintvl,M,cerror,cov);

      if (logps[s] > logpmax)
        { smax = s;
          logpmax = logps[s];
        }
    }

  /* TODO: refactor to something like this:
  for (s = ERROR; s <= DIPLO; s++)
    if ( (logps[s] = calc_logp[s](idx,rintvl,plen,M,perror,cerror,cov)) > logpmax)
      { smax = s;
        logpmax = logps[s];
      }
  */

  bool changed = (rintvl[idx].asgn != smax);

#ifdef DEBUG_REL
  fprintf(stderr,"  logP: ");
  for (s = ERROR; s <= DIPLO; s++)
    fprintf(stderr,"%c=%4.lf%s",stoc[(unsigned char)s],logps[s],(s<DIPLO)?", ":"");
  if (changed)
    fprintf(stderr," *** %c -> %c",stoc[(unsigned char)rintvl[idx].asgn],stoc[(unsigned char)smax]);
  fprintf(stderr,"\n");
#endif

  rintvl[idx].asgn = smax;
  return changed;
}

typedef struct
  { int idx;
    int cnt;
  } Intvl_IC;

static int compare_iic(const void *a, const void *b)
{ return ((Intvl_IC*)a)->cnt - ((Intvl_IC*)b)->cnt;
}

static double logp_e(int idx, Rel_Intvl *rintvl, int plen, P_Error *perror, int cov[])
{ Rel_Intvl ri = rintvl[idx];
  double logp_l, logp_r;
  double logp_po, logp_er;
  
// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"  [E]");
// #endif

  logp_po = logp_poisson(ri.ci,cov[ERROR]);
  logp_er = -INFINITY;
  if (ri.i > 0)
    logp_er = log(perror[ri.i][SELF][DROP]);
  logp_l = MAX(logp_po,logp_er);

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr," [L] logp(PO) = %.1lf, logp(ER) = %.1lf\n",logp_po,logp_er);
// #endif

  logp_po = logp_poisson(ri.cj,cov[ERROR]);
  logp_er = -INFINITY;
  if (ri.j < plen)
    logp_er = log(perror[ri.j][SELF][GAIN]);
  logp_r = MAX(logp_po,logp_er);

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"      [R] logp(PO) = %.1lf, logp(ER) = %.1lf\n",logp_po,logp_er);
// #endif

  return logp_l+logp_r;
}

static double logp_r(int idx, Rel_Intvl *rintvl, int M, int cov[], int pc, int pe)
{ Rel_Intvl ri = rintvl[idx];
  double logp_l, logp_r;
  double logp_sf, logp_er;

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"  [R]");
// #endif

  if (MAX(ri.ci,ri.cj) >= cov[REPEAT])
    { 
// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//       fprintf(stderr," (Larger than R-cov)\n          logp(R) = 0.\n");
// #endif
      return 0.;
    }

  double dr_ratio = 1.+N_SIGMA_R*(1./sqrt(cov[DIPLO]));
  pc = (int)((double)pc*dr_ratio);
  //nc = (int)((double)nc*dr_ratio);

  double _lambda;
  _lambda = (double)pc*(ri.i-pe)/READ_LEN;
  logp_sf = logp_skellam(ri.ci-pc,_lambda);
  logp_er = (pc > ri.ci) ? logp_binom(ri.ci,pc,1-0.01) : -INFINITY;   // TODO: binom test? use ctx
  logp_l = MAX(logp_sf,logp_er);

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr," [L] logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
// #endif

  return logp_l;
}

static double logp_hd(int s, int idx, Rel_Intvl *rintvl, int M, P_Error *cerror, int cov[], int pc, int pe)
{ Rel_Intvl ri = rintvl[idx];
  double logp_l = -INFINITY, logp_r = -INFINITY;
  double logp_sf, logp_er;

  double _lambda = (double)cov[s]*(ri.i-pe)/READ_LEN;
  logp_sf = logp_skellam(ri.ci-pc,_lambda);
  if (pe+1 == ri.i)
    logp_er = log(MAX(cerror[ri.i][OTHERS][DROP],cerror[ri.i][OTHERS][GAIN]));   // TODO: check this prob is what we expect
  else
    logp_er = -INFINITY;
  logp_l = MAX(logp_sf,logp_er);

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr," [L] logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
// #endif

  return logp_l;
}

void classify_reliable(Rel_Intvl *rintvl, int M, Intvl *intvl, int N, int plen, P_Error *perror, P_Error *cerror, int hcov, int dcov)
{ int cov[N_STATE] = {1, dcov+6*(int)sqrt(dcov), hcov, dcov};

  { double dp[M][N_STATE];
    int st[M][N_STATE][2][2];   // H/D, cov/pos
    int bt[M][N_STATE];

    for (int i = 1; i < M; i++)
      for (int s = ERROR; s <= DIPLO; s++)
        { for (int j = 0; j < 2; j++)
            st[i][s][j][0] = -1;
          bt[i][s] = -1;
        }

    for (int s = ERROR; s <= DIPLO; s++)
      { st[0][s][0][0] = hcov;
        st[0][s][0][1] = -1000;
        st[0][s][1][0] = dcov;
        st[0][s][1][1] = -1000;
      }
    dp[0][ERROR] = logp_e(0,rintvl,plen,perror,cov);
    dp[0][REPEAT] = logp_r(0,rintvl,M,cov,dcov,-1000);
    dp[0][HAPLO] = logp_poisson(rintvl[0].ci, hcov);
    st[0][HAPLO][0][0] = rintvl[0].cj;
    st[0][HAPLO][0][1] = rintvl[0].j-1;
    dp[0][DIPLO] = logp_poisson(rintvl[0].ci, dcov);
    st[0][DIPLO][1][0] = rintvl[0].cj;
    st[0][DIPLO][1][1] = rintvl[0].j-1;

#ifdef DEBUG_DP
    fprintf(stderr,"\ndp[0]: ");
    for (int s = ERROR; s <= DIPLO; s++)
      fprintf(stderr,"%c=%.lf, ",stoc[(unsigned char)s],dp[0][s]);
    fprintf(stderr,"\n");
    fprintf(stderr,"st[0]: ");
    for (int s = ERROR; s <= DIPLO; s++)
      fprintf(stderr,"%c=[(%d, %d), (%d, %d)], ",stoc[(unsigned char)s],st[0][s][0][0],st[0][s][0][1],st[0][s][1][0],st[0][s][1][1]);
    fprintf(stderr,"\n");
#endif

    double max_logp, logp;
    int max_s;

    for (int i = 1; i < M; i++)
      { max_logp = -INFINITY;
        max_s = -1;
        for (int s = ERROR; s <= DIPLO; s++)
          { logp = dp[i-1][s];
            if (max_logp < logp)
              { max_logp = logp;
                max_s = s;
              }
          }
        dp[i][ERROR] = max_logp + logp_e(i,rintvl,plen,perror,cov);
        for (int _a = 0; _a < 2; _a++)
          for (int _b = 0; _b < 2; _b++)
            st[i][ERROR][_a][_b] = st[i-1][max_s][_a][_b];
        bt[i][ERROR] = max_s;

        max_logp = -INFINITY;
        max_s = -1;
        for (int s = ERROR; s <= DIPLO; s++)
          { if (st[i-1][s][1][0] == -1)
              continue;
            logp = dp[i-1][s] + logp_r(i,rintvl,M,cov,st[i-1][s][1][0],st[i-1][s][1][1]);
            if (max_logp < logp)
              { max_logp = logp;
                max_s = s;
              }
          }
        dp[i][REPEAT] = max_logp;
        for (int _a = 0; _a < 2; _a++)
          for (int _b = 0; _b < 2; _b++)
            st[i][REPEAT][_a][_b] = st[i-1][max_s][_a][_b];
        bt[i][REPEAT] = max_s;

        max_logp = -INFINITY;
        max_s = -1;
        for (int s = ERROR; s <= DIPLO; s++)
          { if (st[i-1][s][0][0] == -1)
              continue;
            logp = dp[i-1][s] + logp_hd(HAPLO,i,rintvl,M,cerror,cov,st[i-1][s][0][0],st[i-1][s][0][1]);
            if (max_logp < logp)
              { max_logp = logp;
                max_s = s;
              }
          }
        dp[i][HAPLO] = max_logp;
        if (max_s != -1)
          { st[i][HAPLO][0][0] = rintvl[i].cj;
            st[i][HAPLO][0][1] = rintvl[i].j-1;
            for (int _b = 0; _b < 2; _b++)
              st[i][HAPLO][1][_b] = st[i-1][max_s][1][_b];
            bt[i][HAPLO] = max_s;
          }

        max_logp = -INFINITY;
        max_s = -1;
        for (int s = ERROR; s <= DIPLO; s++)
          { if (st[i-1][s][1][0] == -1)
              continue;
            logp = dp[i-1][s] + logp_hd(DIPLO,i,rintvl,M,cerror,cov,st[i-1][s][1][0],st[i-1][s][1][1]);
            if (max_logp < logp)
              { max_logp = logp;
                max_s = s;
              }
          }
        dp[i][DIPLO] = max_logp;
        if (max_s != -1)
          { st[i][DIPLO][1][0] = rintvl[i].cj;
            st[i][DIPLO][1][1] = rintvl[i].j-1;
            for (int _b = 0; _b < 2; _b++)
              st[i][DIPLO][0][_b] = st[i-1][max_s][0][_b];
            bt[i][DIPLO] = max_s;
          }

#ifdef DEBUG_DP
        fprintf(stderr,"\ndp[%d]: ",i);
        for (int s = ERROR; s <= DIPLO; s++)
          fprintf(stderr,"%c=%.lf, ",stoc[(unsigned char)s],dp[i][s]);
        fprintf(stderr,"\n");
        fprintf(stderr,"st[%d]: ",i);
        for (int s = ERROR; s <= DIPLO; s++)
          fprintf(stderr,"%c=[(%d, %d), (%d, %d)], ",stoc[(unsigned char)s],st[i][s][0][0],st[i][s][0][1],st[i][s][1][0],st[i][s][1][1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"bt[%d]: ",i);
        for (int s = ERROR; s <= DIPLO; s++)
          fprintf(stderr,"%c->%c, ",stoc[(unsigned char)bt[i][s]],stoc[(unsigned char)s]);
        fprintf(stderr,"\n");
#endif
      }

    max_logp = -INFINITY;
    max_s = -1;
    for (int s = ERROR; s <= DIPLO; s++)
      { if (max_logp < dp[M-1][s])
          { max_logp = dp[M-1][s];
            max_s = s;
          }
      }
    rintvl[M-1].asgn = max_s;
    for (int i = M-1; i >= 1; i--)
      rintvl[i-1].asgn = bt[i][rintvl[i].asgn];
  }
  
  /*
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
    }*/

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
  bool changed;
  for (counter = 0; counter < REL_MAX_NITER; counter++)
    { changed = false;
      for (int i = M - 1; i >= 0; i--)
        if (update_state(iord[i].idx,rintvl,M,plen,perror,cerror,cov,counter+1))
          changed = true;
      if (!changed) break;

      changed = false;
      for (int i = 0; i < M; i++)
        if (update_state(iord[i].idx,rintvl,M,plen,perror,cerror,cov,-counter-1))
          changed = true;
      if (!changed) break;
    }

#ifdef DEBUG_ITER
  if (counter == REL_MAX_NITER)
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

  // if (pc <= profile[I.i] || profile[I.j-1] >= nc)
  //   return 0.;

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
        logp_er = log(MAX(perror[I.i][OTHERS][DROP],perror[I.i][OTHERS][GAIN]));   // TODO: MIN? MAX?
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
        logp_er = log(MAX(perror[I.j][OTHERS][DROP],perror[I.j][OTHERS][GAIN]));
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
