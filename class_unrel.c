#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "ClassPro.h"

static const int N_SIGMA_R_U           = 3;
static const int N_BASE_EST           = 1000;
static const int N_BASE_EST_MIN           = 1;
static const int N_INTVL_EST           = 1;
static const int N_INTVL_EST_MIN           = 1;
static const double MIN_P_NORMAL           = 0.;

typedef struct
  { int idx;
    int cnt;
  } Intvl_IC;

static int compare_iic(const void *a, const void *b)
{ return ((Intvl_IC*)a)->cnt - ((Intvl_IC*)b)->cnt;
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
