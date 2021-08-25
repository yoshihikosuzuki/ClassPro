#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "ClassPro.h"

static const int N_SIGMA_R_U     = 3;
static const int N_SIGMA_HD_U    = 1;
static double    DR_RATIO_U;
static const int N_BASE_EST      = 1000;
static const int N_BASE_EST_MIN  = 1;
static const int N_INTVL_EST     = 5;
static const int N_INTVL_EST_MIN = 1;

static inline void nn_intvl_u(int idx, Intvl *intvl, int N, char s, int ret[2])
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

static inline void est_cnt_base(int idx, Intvl *intvl, int N, char s, uint16* profile, int ret[2])
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

static inline void est_cnt_intvl_u(int idx, Intvl *intvl, int N, uint16 *profile, char s, int ret[2])
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

static inline double calc_logp_e_u(int idx, Intvl *intvl, uint16 *profile, int plen, P_Error *perror, int cov[])
{ Intvl I = intvl[idx];
  double logp_po, logp_er;

  logp_po = (logp_poisson(profile[I.i],cov[ERROR])
             + logp_poisson(profile[I.j-1],cov[ERROR]));
  // logp_er = MAX(((I.i > 0) ? log(perror[I.i][SELF][DROP]) : -INFINITY),
  //               ((I.j < plen) ? log(perror[I.j][SELF][GAIN]) : -INFINITY));   // TODO: use check_drop/gain?
  logp_er = MAX(I.logpe_i,I.logpe_j);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [E] logp(PO)=%lf, logp(ER)=%lf\n",logp_po,logp_er);
#endif

  return MAX(logp_po,logp_er);
}

static inline double calc_logp_r_u(int idx, Intvl *intvl, int N, uint16 *profile, int plen, int cov[])
{ Intvl I = intvl[idx];
  double logp_er;

  if (MAX(profile[I.i],profile[I.j-1]) >= cov[REPEAT])
    { 
#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
      fprintf(stderr," [R] Larger than R-cov\n");
#endif
      return -10.;
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
  pc = (int)(DR_RATIO_U*pc);
  nc = (int)(DR_RATIO_U*nc);

  if (pc <= profile[I.i] || profile[I.j-1] >= nc)
    return -10.;

  logp_er = (logp_binom(profile[I.i],pc,1-PE_MEAN)   // TODO: binom test? use ctx
             + logp_binom(profile[I.j-1],nc,1-PE_MEAN));

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [R] logp(ER)=%lf\n",logp_er);
#endif

  return logp_er;
}

static inline double calc_logp_hd_u(int s, int idx, Intvl *intvl, int N, uint16 *profile, P_Error *perror, int cov[])
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
        { //fprintf(stderr,"intvl (%d,%d).b: pOD=%lf, pOG=%lf\n",I.i,I.j,perror[I.i][OTHERS][DROP],perror[I.i][OTHERS][GAIN]);
          logp_er = log(MIN(perror[I.i][OTHERS][DROP],perror[I.i][OTHERS][GAIN]));   // TODO: MIN? MAX?
        }
      else
        logp_er = -INFINITY;
      logp_l = MAX(logp_sf,logp_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [%c L] logp(SF)=%lf, logp(ER)=%lf\n",stoc[s],logp_sf,logp_er);
#endif
    }

  if (n < N)
    { double _lambda = (double)cov[s]*(intvl[n].i-I.j+1)/READ_LEN;
      logp_sf = logp_skellam(profile[intvl[n].i]-profile[I.j-1],_lambda);
      if (I.j == intvl[n].i)
        { //fprintf(stderr,"intvl (%d,%d).e: pOD=%lf, pOG=%lf\n",I.i,I.j,perror[I.j][OTHERS][DROP],perror[I.j][OTHERS][GAIN]);
          logp_er = log(MIN(perror[I.j][OTHERS][DROP],perror[I.j][OTHERS][GAIN]));
        }
      else
        logp_er = -INFINITY;
      logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [%c R] logp(SF)=%lf, logp(ER)=%lf\n",stoc[s],logp_sf,logp_er);
#endif
    }

  if (p < 0 || n >= N)
    { if (p < 0 && n >= N)   // TODO: global cov can be wrong
        { logp_l = logp_poisson(profile[I.i],cov[s]);
          logp_r = logp_poisson(profile[I.j-1],cov[s]);
        }
      else if (p < 0)
        logp_l = logp_r;
      else if (n >= N)
        logp_r = logp_l;
    }

  return logp_l+logp_r;
}

static inline double calc_logp_h_u(int idx, Intvl *intvl, int N, uint16 *profile, P_Error *perror, int cov[])
{ Intvl I = intvl[idx];
  int ibc = profile[I.i], iec = profile[I.j-1];
  
  // H must not exceed D - n_sigma
  int nn_idx[2];
  nn_intvl_u(idx,intvl,N,'D',nn_idx);
  int p = nn_idx[0];
  int n = nn_idx[1];
  if (p >= 0)
    { int pc = profile[intvl[p].j-1];
      if (minus_sigma(pc,N_SIGMA_HD_U) < ibc)
        return -INFINITY;
    }
  if (n < N)
    { int nc = profile[intvl[n].i];
      if (iec > minus_sigma(nc,N_SIGMA_HD_U))
        return -INFINITY;
    }

  return calc_logp_hd_u(HAPLO,idx,intvl,N,profile,perror,cov);
}

static inline double calc_logp_d_u(int idx, Intvl *intvl, int N, uint16 *profile, P_Error *perror, int cov[])
{ return calc_logp_hd_u(DIPLO,idx,intvl,N,profile,perror,cov);
}

static void update_state_u(int idx, int d, Intvl *intvl, int N, uint16 *profile, int plen, P_Error *perror, int cov[])
{ double logp, logpmax = -INFINITY;
  int smax = -1;

#ifdef DEBUG_UNREL
      fprintf(stderr,"\n<%d> Intvl[%d] @ (%d,%d) %d -> %d: class=%c\n",d,idx,intvl[idx].i,intvl[idx].j,profile[intvl[idx].i],profile[intvl[idx].j-1],(intvl[idx].asgn == N_STATE) ? '-' : stoc[(int)intvl[idx].asgn]);
#endif

  for (int s = ERROR; s <= DIPLO; s++)
    { if (s == ERROR)
        logp = calc_logp_e_u(idx,intvl,profile,plen,perror,cov);
      else if (s == REPEAT)
        logp = calc_logp_r_u(idx,intvl,N,profile,plen,cov);
      else if (s == HAPLO)
        logp = calc_logp_h_u(idx,intvl,N,profile,perror,cov);
      else
        logp = calc_logp_d_u(idx,intvl,N,profile,perror,cov);

      if (logp > logpmax)
        { logpmax = logp;
          smax = s;
        }

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
      fprintf(stderr,"logp(%c)=%lf\n",stoc[s],logp);
#endif
    }

#ifdef DEBUG
  if (smax == -1)
    { fprintf(stderr,"No valid probability for interval %d\n",idx);
      exit(1);
    }
#endif

  if (intvl[idx].asgn != smax)
    {
#ifdef DEBUG_UNREL
      fprintf(stderr,"state updated @ %d: %c -> %c\n",idx,(intvl[idx].asgn == N_STATE) ? '-' : stoc[(int)intvl[idx].asgn],stoc[smax]);
#endif

      intvl[idx].asgn = smax;
    }

  return;
}

typedef struct
  { int idx;
    int cnt;
  } Intvl_IC;

static int compare_iic(const void *a, const void *b)
{ return ((Intvl_IC*)a)->cnt - ((Intvl_IC*)b)->cnt;
}

// void extend_e()
// { 
// }

void classify_unreliable(uint16 *profile, int plen, Intvl *intvl, int N, P_Error *perror, int hcov, int dcov)
{ int cov[N_STATE] = {1, dcov+6*(int)sqrt(dcov), hcov, dcov};
  DR_RATIO_U = 1.+N_SIGMA_R_U*(1./sqrt(cov[DIPLO]));

  // If no H/D-mers, then do initial classification
  // int local_cov[N_STATE] = {1, }

  bool is_fixed[N];
  for (int i = 0; i < N; i++)
    is_fixed[i] = (intvl[i].is_rel && (intvl[i].asgn == HAPLO || intvl[i].asgn == DIPLO));

  Intvl_IC iord[N];
  for (int i = 0; i < N; i++)
    { iord[i].idx = i;
      iord[i].cnt = MIN(profile[intvl[i].i],profile[intvl[i].j-1]);
    }
  qsort(iord,N,sizeof(Intvl_IC),compare_iic);

  for (int i = N-1; i >= 0; i--)
    if (!is_fixed[iord[i].idx])
      update_state_u(iord[i].idx,0,intvl,N,profile,plen,perror,cov);

  for (int i = 0; i < N; i++)
    if (!is_fixed[iord[i].idx])
      update_state_u(iord[i].idx,1,intvl,N,profile,plen,perror,cov);

#ifdef DEBUG_ITER
  fprintf(stderr,"         ");
  for (int i = 0; i < N; i++)
    { if (!is_fixed[i])
        fprintf(stderr,"+");
      else
        fprintf(stderr," ");
    }
  fprintf(stderr,"\n  Intvl: ");
  for (int i = 0; i < N; i++)
    fprintf(stderr,"%c",stoc[(int)intvl[i].asgn]);
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

#if defined(DEBUG_ITER) && defined(DEBUG_SINGLE)
  char asgn[plen];
  for (int i = 0; i < N; i++)
    for (int j = intvl[i].i; j < intvl[i].j; j++)
      asgn[j] = intvl[i].asgn;

  fprintf(stderr,"  Final: ");
  for (int i = 0; i < plen; i++)
    fprintf(stderr,"%c",stoc[(int)asgn[i]]);
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
