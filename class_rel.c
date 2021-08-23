#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "ClassPro.h"

static const int    N_SIGMA_R = 3;
static const int    N_SIGMA_FLUC = 1;
static double       DR_RATIO_R;
static const double THRES_LOGP_FLUC = log(0.001);   // TODO: make `p_skellam` and use p?
static const double THRES_LOGP_REP = -50.;
static const double LOGP_E_BASE = -10.;
static const double LOGP_R = -10.;
static const double PE_MEAN = 0.01;
static const int    POS_OFFSET = 1000;

/*******************************************************************************************
 *
 *  Utilities
 *
 ********************************************************************************************/

static inline int plus_sigma(int cnt, int n_sigma)
{ return cnt + (int)(n_sigma * sqrt(cnt));
}

static inline int minus_sigma(int cnt, int n_sigma)
{ return cnt - (int)(n_sigma * sqrt(cnt));
}

static inline double logp_fluctuation(int i, int j, int ci, int cj, int mean_cov)
{ double _lambda = (double)mean_cov*(j-i)/READ_LEN;
  double logp = logp_skellam(cj-ci,_lambda);
  return logp;
}

static inline int find_nn_fw(int i, int s, char *asgn, int end)
{ int idx = i;
  while (idx < end && asgn[idx] != s)
    idx++;
  return idx;
}

static inline int find_nn_bw(int i, int s, char *asgn, int start)
{ int idx = i;
  while (idx >= start && asgn[idx] != s)
    idx--;
  return idx;
}

static inline double calc_d_to_h_fw(int i, char *asgn, Rel_Intvl *rintvl)
{ // Find D2 -> H -> D1 within [0..i]
  int d1_idx, h_idx, d2_idx;
  int d1_pos, h_pos, d2_pos;
  int d1_cnt, h_cnt, d2_cnt;
  double pos_ratio, hd_ratio;

  if ((d1_idx = find_nn_bw(i,DIPLO,asgn,0)) < 0)
    return -INFINITY;
  if ((h_idx = find_nn_bw(d1_idx-1,HAPLO,asgn,0)) < 0)
    return -INFINITY;
  if ((d2_idx = find_nn_bw(h_idx-1,DIPLO,asgn,0)) < 0)
    return -INFINITY;
  d1_idx = find_nn_fw(h_idx+1,DIPLO,asgn,i+1);

  d2_pos = rintvl[d2_idx].j-1;   // TODO: mid point is better?
  d2_cnt = rintvl[d2_idx].cj;
  h_pos = rintvl[h_idx].j-1;
  h_cnt = rintvl[h_idx].cj;
  d1_pos = rintvl[d1_idx].i;
  d1_cnt = rintvl[d1_idx].ci;

  pos_ratio = (double)(h_pos-d2_pos)/(d1_pos-d2_pos);
  hd_ratio = (pos_ratio*(d1_cnt-d2_cnt)+d2_cnt)/h_cnt;

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"D[%d].e=%d@%d -> H[%d].b=%d@%d -> D[%d].b=%d@%d pos_ratio=%lf hd_ratio=%lf\n",d2_idx,d2_cnt,d2_pos,h_idx,h_cnt,h_pos,d1_idx,d1_cnt,d1_pos,pos_ratio,hd_ratio);
// #endif

  return hd_ratio;
}

static inline double calc_h_to_d_fw(int i, char *asgn, Rel_Intvl *rintvl)
{ // Find H2 -> D -> H1 within [0..i]
  int h1_idx, d_idx, h2_idx;
  int h1_pos, d_pos, h2_pos;
  int h1_cnt, d_cnt, h2_cnt;
  double pos_ratio, hd_ratio;

  if ((h1_idx = find_nn_bw(i,HAPLO,asgn,0)) < 0)
    return -INFINITY;
  if ((d_idx = find_nn_bw(h1_idx-1,DIPLO,asgn,0)) < 0)
    return -INFINITY;
  if ((h2_idx = find_nn_bw(d_idx-1,HAPLO,asgn,0)) < 0)
    return -INFINITY;
  h1_idx = find_nn_fw(d_idx+1,HAPLO,asgn,i+1);

  h2_pos = rintvl[h2_idx].j-1;   // TODO: mid point is better?
  h2_cnt = rintvl[h2_idx].cj;
  d_pos = rintvl[d_idx].j-1;
  d_cnt = rintvl[d_idx].cj;
  h1_pos = rintvl[h1_idx].i;
  h1_cnt = rintvl[h1_idx].ci;

  pos_ratio = (double)(d_pos-h2_pos)/(h1_pos-h2_pos);
  hd_ratio = (double)d_cnt/(pos_ratio*(h1_cnt-h2_cnt)+h2_cnt);

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"H[%d].e=%d@%d -> D[%d].e=%d@%d -> H[%d].b=%d@%d pos_ratio=%lf hd_ratio=%lf\n",h2_idx,h2_cnt,h2_pos,d_idx,d_cnt,d_pos,h1_idx,h1_cnt,h1_pos,pos_ratio,hd_ratio);
// #endif

  return hd_ratio;
}

static inline double calc_d_to_h_bw(int i, char *asgn, Rel_Intvl *rintvl, int M)
{ // Find D2 -> H -> D1 within [i..M-1]
  int d1_idx, h_idx, d2_idx;
  int d1_pos, h_pos, d2_pos;
  int d1_cnt, h_cnt, d2_cnt;
  double pos_ratio, hd_ratio;

  if ((d1_idx = find_nn_fw(i,DIPLO,asgn,M)) >= M)
    return -INFINITY;
  if ((h_idx = find_nn_fw(d1_idx+1,HAPLO,asgn,M)) >= M)
    return -INFINITY;
  if ((d2_idx = find_nn_fw(h_idx+1,DIPLO,asgn,M)) >= M)
    return -INFINITY;
  d1_idx = find_nn_bw(h_idx-1,DIPLO,asgn,i);

  d2_pos = rintvl[d2_idx].i;   // TODO: mid point is better?
  d2_cnt = rintvl[d2_idx].ci;
  h_pos = rintvl[h_idx].i;
  h_cnt = rintvl[h_idx].ci;
  d1_pos = rintvl[d1_idx].j-1;
  d1_cnt = rintvl[d1_idx].cj;

  pos_ratio = (double)(h_pos-d1_pos)/(d2_pos-d1_pos);
  hd_ratio = (pos_ratio*(d2_cnt-d1_cnt)+d1_cnt)/h_cnt;

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"D[%d].e=%d@%d <- H[%d].b=%d@%d <- D[%d].b=%d@%d pos_ratio=%lf hd_ratio=%lf\n",d1_idx,d1_cnt,d1_pos,h_idx,h_cnt,h_pos,d2_idx,d2_cnt,d2_pos,pos_ratio,hd_ratio);
// #endif

  return hd_ratio;
}

static inline double calc_h_to_d_bw(int i, char *asgn, Rel_Intvl *rintvl, int M)
{ // Find H2 -> D -> H1 within [0..i]
  int h1_idx, d_idx, h2_idx;
  int h1_pos, d_pos, h2_pos;
  int h1_cnt, d_cnt, h2_cnt;
  double pos_ratio, hd_ratio;

  if ((h1_idx = find_nn_fw(i,HAPLO,asgn,M)) >= M)
    return -INFINITY;
  if ((d_idx = find_nn_fw(h1_idx+1,DIPLO,asgn,M)) >= M)
    return -INFINITY;
  if ((h2_idx = find_nn_fw(d_idx+1,HAPLO,asgn,M)) >= M)
    return -INFINITY;
  h1_idx = find_nn_bw(d_idx-1,HAPLO,asgn,i);

  h2_pos = rintvl[h2_idx].i;   // TODO: mid point is better?
  h2_cnt = rintvl[h2_idx].ci;
  d_pos = rintvl[d_idx].i;
  d_cnt = rintvl[d_idx].ci;
  h1_pos = rintvl[h1_idx].j-1;
  h1_cnt = rintvl[h1_idx].cj;

  pos_ratio = (double)(d_pos-h1_pos)/(h2_pos-h1_pos);
  hd_ratio = (double)d_cnt/(pos_ratio*(h2_cnt-h1_cnt)+h1_cnt);

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"H[%d].e=%d@%d <- D[%d].e=%d@%d <- H[%d].b=%d@%d pos_ratio=%lf hd_ratio=%lf\n",h1_idx,h1_cnt,h1_pos,d_idx,d_cnt,d_pos,h2_idx,h2_cnt,h2_pos,pos_ratio,hd_ratio);
// #endif

  return hd_ratio;
}

/*******************************************************************************************
 *
 *  Transition probabilities
 *
 ********************************************************************************************/

static inline double logp_e_fw(int idx, Rel_Intvl *rintvl, int plen, P_Error *perror, int cov[])
{ Rel_Intvl ri = rintvl[idx];
  double logp_l, logp_r;
  double logp_po, logp_er;

  logp_po = logp_poisson(ri.ci,cov[ERROR])+LOGP_E_BASE;
  // logp_po = -INFINITY;
  logp_er = -INFINITY;
  if (ri.i > 0)
    logp_er = log(perror[ri.i][SELF][DROP]);
  logp_l = MAX(logp_po,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [E L] logp(PO) = %.1lf, logp(ER) = %.1lf\n",logp_po,logp_er);
#endif

  logp_po = logp_poisson(ri.cj,cov[ERROR])+LOGP_E_BASE;
  // logp_po = -INFINITY;
  logp_er = -INFINITY;
  if (ri.j < plen)
    logp_er = log(perror[ri.j][SELF][GAIN]);
  logp_r = MAX(logp_po,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [E R] logp(PO) = %.1lf, logp(ER) = %.1lf\n",logp_po,logp_er);
#endif

  return logp_l+logp_r;
}

static inline double logp_e_bw(int idx, Rel_Intvl *rintvl, int plen, P_Error *perror, int cov[])
{ return logp_e_fw(idx,rintvl,plen,perror,cov);
}

static inline double logp_r_fw(int idx, Rel_Intvl *rintvl, int cov[], int st[])
{ Rel_Intvl ri = rintvl[idx];
  int pc = st[0], pe = st[1];

  double logp_l;
  double logp_sf, logp_er;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [R L] ");
#endif

  if (MAX(ri.ci,ri.cj) >= cov[REPEAT])
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,"(Larger than R-cov)\n          logp(R) = 0.\n");
#endif
      return LOGP_R;
    }

  if (MAX(ri.ci,ri.cj) >= pc)
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,"(Larger than imaginary R-cov = %d)\n          logp(R) = 0.\n",pc);
#endif
      return LOGP_R;
    }

  double _lambda = (double)pc*(ri.i-pe)/READ_LEN;
  logp_sf = logp_skellam(ri.ci-pc,_lambda);
  logp_er = (pc > ri.ci) ? logp_binom(ri.ci,pc,1-PE_MEAN) : -INFINITY;   // TODO: binom test? use ctx
  // logp_er = logp_poisson(ri.ci,pc)
  logp_l = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
#endif

  return logp_l;
}

static inline double logp_r_bw(int idx, Rel_Intvl *rintvl, int cov[], int st[])
{ Rel_Intvl ri = rintvl[idx];
  int nc = st[0], nb = st[1];

  double logp_r;
  double logp_sf, logp_er;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [R R] ");
#endif

  if (MAX(ri.ci,ri.cj) >= cov[REPEAT])
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr," (Larger than R-cov)\n          logp(R) = 0.\n");
#endif
      return LOGP_R;
    }

  if (MAX(ri.ci,ri.cj) >= nc)
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,"(Larger than imaginary R-cov = %d)\n          logp(R) = 0.\n",nc);
#endif
      return LOGP_R;
    }

  double _lambda = (double)nc*(nb-ri.j+1)/READ_LEN;
  logp_sf = logp_skellam(nc-ri.cj,_lambda);
  logp_er = (nc > ri.cj) ? logp_binom(ri.cj,nc,1-PE_MEAN) : -INFINITY;   // TODO: binom test? use ctx
  logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
#endif

  return logp_r;
}

static inline double logp_hd_fw(int s, int idx, Rel_Intvl *rintvl, P_Error *cerror, int cov[], int st[])
{ Rel_Intvl ri = rintvl[idx];
  int pc = st[0], pe = st[1];

  double logp_l;
  double logp_sf, logp_er;

  double _lambda = (double)cov[s]*(ri.i-pe)/READ_LEN;
  logp_sf = logp_skellam(ri.ci-pc,_lambda);
  logp_er = -INFINITY;   // TODO: do not use count at drop/gain by errors in others as H/D-cov
  // if (pe+1 == ri.i)
  //   logp_er = log(MAX(cerror[ri.i][OTHERS][DROP],cerror[ri.i][OTHERS][GAIN]));   // TODO: check this prob is what we expect
  logp_l = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr," [%c L] logp(SF) = %.1lf, logp(ER) = %.1lf\n",stoc[(unsigned char)s],logp_sf,logp_er);
#endif

  return logp_l;
}

static inline double logp_hd_bw(int s, int idx, Rel_Intvl *rintvl, P_Error *cerror, int cov[], int st[])
{ Rel_Intvl ri = rintvl[idx];
  int nc = st[0], nb = st[1];

  double logp_r;
  double logp_sf, logp_er;

  double _lambda = (double)cov[s]*(nb-ri.j+1)/READ_LEN;
  logp_sf = logp_skellam(nc-ri.cj,_lambda);
  logp_er = -INFINITY;   // TODO: do not use count at drop/gain by errors in others as H/D-cov
  // if (ri.j == nb)
  //   logp_er = log(MAX(cerror[ri.j][OTHERS][DROP],cerror[ri.j][OTHERS][GAIN]));   // TODO: check this prob is what we expect
  logp_r = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr," [%c R] logp(SF) = %.1lf, logp(ER) = %.1lf\n",stoc[(unsigned char)s],logp_sf,logp_er);
#endif

  return logp_r;
}

static inline double logp_h_fw(int idx, Rel_Intvl *rintvl, P_Error *cerror, int cov[], int st[])
{ return logp_hd_fw(HAPLO,idx,rintvl,cerror,cov,st);
}

static inline double logp_h_bw(int idx, Rel_Intvl *rintvl, P_Error *cerror, int cov[], int st[])
{ return logp_hd_bw(HAPLO,idx,rintvl,cerror,cov,st);
}

static inline double logp_d_fw(int idx, Rel_Intvl *rintvl, P_Error *cerror, int cov[], int st[])
{ return logp_hd_fw(DIPLO,idx,rintvl,cerror,cov,st);
}

static inline double logp_d_bw(int idx, Rel_Intvl *rintvl, P_Error *cerror, int cov[], int st[])
{ return logp_hd_bw(DIPLO,idx,rintvl,cerror,cov,st);
}

/*******************************************************************************************
 *
 *  Forward DP
 *
 ********************************************************************************************/

static void forward_update(int i, double **dp, int ****st, char ***bt, bool *rep, int POS_INIT, Rel_Intvl *rintvl, int plen, P_Error *perror, P_Error *cerror, int cov[])
{ Rel_Intvl ri = rintvl[i];

  double logp_trans[N_STATE][N_STATE];   // s @ i-1 -> t @ i
  double logp, max_logp, psum, logp_fluc;
  int max_s;
  int **prev_st;
  int pc, pe, icms, icps;
  bool hr_violated, hd_violated, dr_violated;

  for (int s = ERROR; s <= DIPLO; s++)
    for (int t = ERROR; t <= DIPLO; t++)
      logp_trans[s][t] = -INFINITY;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"\n### [Intvl %d] @ (%d, %d): %d -> %d\n",i,ri.i,ri.j,ri.ci,ri.cj);
#endif

  // Compute and normalize transition probabilities so that \sum_{t'} p(s -> t') = 1
  psum = 0.;
  for (int s = ERROR; s <= DIPLO; s++)
    { if (dp[i-1][s] == -INFINITY)
        continue;
      prev_st = st[i-1][s];
      for (int t = ERROR; t <= DIPLO; t++)
        { if (t == ERROR)
            { logp = logp_e_fw(i,rintvl,plen,perror,cov);
            }
          else if (t == REPEAT)
            { // prev H < this R
              icms = minus_sigma(ri.ci,N_SIGMA_FLUC);
              pc = prev_st[HAPLO][0];
              pe = prev_st[HAPLO][1];
              logp_fluc = logp_fluctuation(pe-1,ri.i,pc,icms,cov[HAPLO]);
              hr_violated = !(pe == POS_INIT) && !(pc < icms) && !(logp_fluc > THRES_LOGP_FLUC);

              // prev D < this R
              pc = prev_st[DIPLO][0];
              pe = prev_st[DIPLO][1];
              logp_fluc = logp_fluctuation(pe-1,ri.i,pc,icms,cov[DIPLO]);
              // fprintf(stderr,"pe=%d,POS_INIT=%d,pc=%d,icms=%d, logp_fluc (%d @ %d-> %d @ %d)=%lf, THRES=%lf\n",pe,POS_INIT,pc,icms,pe-1,pc,ri.i,icms,logp_fluc,THRES_LOGP_FLUC);
              dr_violated = !(pe == POS_INIT) && !(pc < icms) && !(logp_fluc > THRES_LOGP_FLUC);

              if (hr_violated || dr_violated)
                logp = -INFINITY;
              else
                logp = logp_r_fw(i,rintvl,cov,prev_st[REPEAT]);
            }
          else if (t == HAPLO)
            { // this H < prev D
              icps = plus_sigma(ri.ci,N_SIGMA_FLUC);
              pc = prev_st[DIPLO][0];
              pe = prev_st[DIPLO][1];
              logp_fluc = logp_fluctuation(pe-1,ri.i,pc,icps,cov[DIPLO]);
              hd_violated = !(pe == POS_INIT) && !(icps < pc) && !(logp_fluc > THRES_LOGP_FLUC);

              if (hd_violated)
                logp = -INFINITY;
              else
                logp = logp_h_fw(i,rintvl,cerror,cov,prev_st[HAPLO]);
            }
          else   // DIPLO
            { // prev H < this D
              icms = minus_sigma(ri.ci,N_SIGMA_FLUC);
              pc = prev_st[HAPLO][0];
              pe = prev_st[HAPLO][1];
              logp_fluc = logp_fluctuation(pe-1,ri.i,pc,icms,cov[HAPLO]);
              hd_violated = !(pe == POS_INIT) && !(pc < icms) && !(logp_fluc > THRES_LOGP_FLUC);

              // this H < prev D
              icps = plus_sigma(ri.ci,N_SIGMA_FLUC);
              pc = prev_st[REPEAT][0];
              pe = prev_st[REPEAT][1];
              logp_fluc = logp_fluctuation(pe-1,ri.i,pc,icps,cov[REPEAT]);
              dr_violated = !(pe == POS_INIT) && !(icps < pc) && !(logp_fluc > THRES_LOGP_FLUC);

              if (hd_violated || dr_violated)
                logp = -INFINITY;
              else
                logp = logp_d_fw(i,rintvl,cerror,cov,prev_st[DIPLO]);
            }

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
          fprintf(stderr,"%c -> %c: logp = %.lf\n",stoc[s],stoc[t],logp);
#endif

          logp_trans[s][t] = exp(logp);
          psum += logp_trans[s][t];
        }
    }
  for (int s = ERROR; s <= DIPLO; s++)
    { if (dp[i-1][s] == -INFINITY)
        continue;
      for (int t = ERROR; t <= DIPLO; t++)
        logp_trans[s][t] = (psum > 0.) ? log(logp_trans[s][t]/psum) : -INFINITY;
    }

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"Normalized logp_trans:\n");
//   for (int s = ERROR; s <= DIPLO; s++)
//     for (int t = ERROR; t <= DIPLO; t++)
//       fprintf(stderr,"  %c -> %c: logp = %.lf\n",stoc[s],stoc[t],logp_trans[s][t]);
// #endif

  bool only_r = true;
  for (int t = ERROR; t <= DIPLO; t++)
    { max_logp = -INFINITY;
      max_s = -1;
      for (int s = ERROR; s <= DIPLO; s++)
        { if (dp[i-1][s] == -INFINITY)
            continue;
          logp = dp[i-1][s]+logp_trans[s][t];
          if (max_logp < logp)
            { max_logp = logp;
              max_s = s;
              if (t != REPEAT && logp_trans[s][t] > THRES_LOGP_REP)
                only_r = false;
            }
        }
    }

  if (only_r)
    { 
#ifdef DEBUG_DP
      fprintf(stderr,"Absolutely REPEAT @ %d\n",i);
#endif
      rep[i] = true;
      for (int s = ERROR; s <= DIPLO; s++)
        { dp[i][s] = dp[i-1][s];
          if (dp[i-1][s] == -INFINITY)
            continue;
          for (int ii = 0; ii < i; ii++)
            bt[i][s][ii] = bt[i-1][s][ii];
          bt[i][s][i] = s;
          for (int t = ERROR; t <= DIPLO; t++)
            for (int n = 0; n < 2; n++)
              st[i][s][t][n] = st[i-1][s][t][n];
        }
      // NOTE: DO NOT USE `ri`
      rintvl[i].i = rintvl[i-1].i;
      rintvl[i].j = rintvl[i-1].j;
      rintvl[i].ci = rintvl[i-1].ci;
      rintvl[i].cj = rintvl[i-1].cj;
      rintvl[i].asgn = rintvl[i-1].asgn;
      return;
    }

  int curr_h, curr_d, curr_r;
  double dh_ratio;

  for (int t = ERROR; t <= DIPLO; t++)
    { max_logp = -INFINITY;
      max_s = -1;
      for (int s = ERROR; s <= DIPLO; s++)
        { if (dp[i-1][s] == -INFINITY)
            continue;
          logp = dp[i-1][s]+logp_trans[s][t];
          if (max_logp < logp)
            { max_logp = logp;
              max_s = s;
            }
        }
      dp[i][t] = max_logp;
      if (max_s != -1)
        { for (int ii = 0; ii < i; ii++)
            bt[i][t][ii] = bt[i-1][max_s][ii];
          bt[i][t][i] = t;

          if (t == ERROR)
            { for (int s = REPEAT; s <= DIPLO; s++)
                for (int n = 0; n < 2; n++)
                  st[i][t][s][n] = st[i-1][max_s][s][n];
            }
          else if (t == REPEAT)
            { for (int s = HAPLO; s <= DIPLO; s++)
                { st[i][t][s][0] = st[i-1][max_s][s][0];
                  st[i][t][s][1] = ri.j-1-POS_OFFSET;
                }
              if (st[i-1][max_s][REPEAT][0] < MIN(cov[REPEAT],ri.cj))
                for (int n = 0; n < 2; n++)
                  st[i][t][REPEAT][n] = st[i-1][max_s][REPEAT][n];
              else
                { st[i][t][REPEAT][0] = MIN(cov[REPEAT],ri.cj);
                  st[i][t][REPEAT][1] = ri.j-1;
                }
            }
          else if (t == HAPLO)
            { curr_h = ri.cj;
              dh_ratio = calc_h_to_d_fw(i,bt[i][t],rintvl);
              // fprintf(stderr,"i=%d,max_s=%d,t=%d,dh_ratio=%lf\n",i,max_s,t,dh_ratio);
              if (dh_ratio == -INFINITY)
                { bool has_d = false;
                  for (int ii = 0; ii < i; ii++)
                    if (bt[i][t][ii] == DIPLO)
                      has_d = true;
                  if (has_d)
                    curr_d = st[i-1][max_s][DIPLO][0];   // TODO: update pos as well?
                  else
                    curr_d = MAX(cov[DIPLO],curr_h*2);
                }
              else
                curr_d = (int)(dh_ratio*curr_h);
              curr_r = (int)(DR_RATIO_R*curr_d);
              st[i][t][HAPLO][0] = curr_h;
              st[i][t][HAPLO][1] = ri.j-1;
              st[i][t][DIPLO][0] = curr_d;
              st[i][t][DIPLO][1] = ri.j-1-POS_OFFSET;
              st[i][t][REPEAT][0] = curr_r;
              st[i][t][REPEAT][1] = ri.j-1-POS_OFFSET;
            }
          else   // DIPLO
            { curr_d = ri.cj;
              dh_ratio = calc_d_to_h_fw(i,bt[i][t],rintvl);
              if (dh_ratio == -INFINITY)
                { bool has_h = false;
                  for (int ii = 0; ii < i; ii++)
                    if (bt[i][t][ii] == HAPLO)
                      has_h = true;
                  if (has_h)
                    curr_h = st[i-1][max_s][HAPLO][0];   // TODO: update pos as well?
                  else
                    curr_h = curr_d/2;
                }
              else
                curr_h = (int)((double)curr_d/dh_ratio);
              curr_r = (int)(DR_RATIO_R*curr_d);
              st[i][t][HAPLO][0] = curr_h;
              st[i][t][HAPLO][1] = ri.j-1-POS_OFFSET;
              st[i][t][DIPLO][0] = curr_d;
              st[i][t][DIPLO][1] = ri.j-1;
              st[i][t][REPEAT][0] = curr_r;
              st[i][t][REPEAT][1] = ri.j-1-POS_OFFSET;
            }
        }

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      for (int s = ERROR; s <= DIPLO; s++)
        { logp = dp[i-1][s]+logp_trans[s][t];
          fprintf(stderr,"(%c ->) %c: logp = %.lf\n",stoc[s],stoc[t],logp);
        }
#endif
    }

  return;
}

static inline void forward_print(int i, double **dp, int ****st, char ***bt, Rel_Intvl *rintvl)
{ Rel_Intvl ri = rintvl[i];
  double max_logp = -INFINITY;
  int max_s = -1;

  for (int s = ERROR; s <= DIPLO; s++)
    { if (max_logp < dp[i][s])
        { max_logp = dp[i][s];
          max_s = s;
        }
    }

  fprintf(stderr,"\n@@@ Intvl[%d] (%d, %d): %d -> %d\n",i,ri.i,ri.j,ri.ci,ri.cj);
  for (int s = ERROR; s <= DIPLO; s++)
    { if (i == 0)
        { fprintf(stderr,"  %c: logp=%5.lf, H/D/R=(%2d,%2d,%3d) %s\n",stoc[s],dp[i][s],st[i][s][HAPLO][0],st[i][s][DIPLO][0],st[i][s][REPEAT][0],(s == max_s) ? "***" : "");
        }
      else
        { int prev_s = bt[i][s][i-1];
          if (prev_s == ' ')   // impossible state
            continue;
          int **prev_st = st[i-1][prev_s];
          fprintf(stderr,"  (%2d,%2d,%3d; %c->) %c: logp=%5.lf, H/D/R=(%2d,%2d,%3d) ",prev_st[HAPLO][0],prev_st[DIPLO][0],prev_st[REPEAT][0],stoc[prev_s],stoc[s],dp[i][s],st[i][s][HAPLO][0],st[i][s][DIPLO][0],st[i][s][REPEAT][0]); 
          for (int ii = 0; ii <= i; ii++)
            fprintf(stderr,"%c",stoc[(int)bt[i][s][ii]]);
          fprintf(stderr," %s\n",(s == max_s) ? "***" : "");
        }
    }
  fprintf(stderr,"\n");
  
  return;
}

static void forward_dp(double **dp, int ****st, char ***bt, bool *rep, char *asgn, Rel_Intvl *_rintvl, int M, int plen, P_Error *perror, P_Error *cerror, int cov[])
{ // Copy rintvl
  Rel_Intvl *rintvl = Malloc(sizeof(Rel_Intvl)*M,"Copy rintvl");
  for (int i = 0; i < M; i++)
    { rintvl[i].i = _rintvl[i].i;
      rintvl[i].j = _rintvl[i].j;
      rintvl[i].ci = _rintvl[i].ci;
      rintvl[i].cj = _rintvl[i].cj;
      rintvl[i].asgn = _rintvl[i].asgn;
    }
  
  // Initialization
  for (int i = 0; i < M; i++)
    { rep[i] = false;
      for (int s = ERROR; s <= DIPLO; s++)
        { dp[i][s] = -INFINITY;
          for (int t = REPEAT; t <= DIPLO; t++)
            st[i][s][t][0] = -1;
          for (int j = 0; j < M; j++)
            bt[i][s][j] = ' ';
        }
    }

  const int POS_INIT = -POS_OFFSET;
  int i;

  // First interval
  i = 0;
  Rel_Intvl ri = rintvl[i];
  for (int s = ERROR; s <= DIPLO; s++)
    { for (int t = REPEAT; t <= DIPLO; t++)
        { st[i][s][t][0] = cov[t];
          st[i][s][t][1] = POS_INIT;
        }
      bt[i][s][i] = s;
    }
  
  dp[i][ERROR] = logp_e_fw(i,rintvl,plen,perror,cov);
  dp[i][REPEAT] = logp_r_fw(i,rintvl,cov,st[i][REPEAT][REPEAT]);
  
  dp[i][HAPLO] = logp_poisson(ri.ci, cov[HAPLO]);
  st[i][HAPLO][HAPLO][0] = ri.cj;
  st[i][HAPLO][HAPLO][1] = ri.j-1;
  st[i][HAPLO][DIPLO][0] = ri.cj*2;
  st[i][HAPLO][DIPLO][1] = ri.j-1-POS_OFFSET;
  st[i][HAPLO][REPEAT][0] = (int)(DR_RATIO_R*2*ri.cj);
  st[i][HAPLO][REPEAT][1] = ri.j-1-POS_OFFSET;
  
  
  dp[i][DIPLO] = logp_poisson(ri.ci, cov[DIPLO]);
  st[i][DIPLO][HAPLO][0] = ri.cj/2;
  st[i][DIPLO][HAPLO][1] = ri.j-1-POS_OFFSET;
  st[i][DIPLO][DIPLO][0] = ri.cj;
  st[i][DIPLO][DIPLO][1] = ri.j-1;
  st[i][DIPLO][REPEAT][0] = (int)(DR_RATIO_R*ri.cj);
  st[i][DIPLO][REPEAT][1] = ri.j-1-POS_OFFSET;

  double psum = 0.;
  for (int s = ERROR; s <= DIPLO; s++)
    psum += exp(dp[i][s]);
  for (int s = ERROR; s <= DIPLO; s++)
    dp[i][s] = log(exp(dp[i][s])/psum);   // TODO: reduce log/exp

#ifdef DEBUG_DP
  forward_print(i,dp,st,bt,rintvl);
#endif

  // The others
  for (i = 1; i < M; i++)
    { forward_update(i,dp,st,bt,rep,POS_INIT,rintvl,plen,perror,cerror,cov);

#ifdef DEBUG_DP
      forward_print(i,dp,st,bt,rintvl);
#endif
    }

  // Traceback
  double max_logp;
  int max_s;

  max_logp = -INFINITY;
  max_s = -1;
  for (int s = ERROR; s <= DIPLO; s++)
    { if (max_logp < dp[M-1][s])
        { max_logp = dp[M-1][s];
          max_s = s;
        }
    }
  for (i = 0; i < M; i++)
    asgn[i] = (rep[i]) ? REPEAT : bt[M-1][max_s][i];

  free(rintvl);
  return;
}

/*******************************************************************************************
 *
 *  Backward DP
 *
 ********************************************************************************************/

static void backward_update(int i, double **dp, int ****st, char ***bt, bool *rep, int POS_INIT, Rel_Intvl *rintvl, int M, int plen, P_Error *perror, P_Error *cerror, int cov[])
{ Rel_Intvl ri = rintvl[i];
  int ic = ri.cj;   // TODO: integrate forward and backward?

  double logp_trans[N_STATE][N_STATE];   // t @ i <- s @ i+1
  double logp, max_logp, psum, logp_fluc;
  int max_s;
  int **next_st;
  int nc, nb, icms, icps;
  bool hr_violated, hd_violated, dr_violated;

  for (int s = ERROR; s <= DIPLO; s++)
    for (int t = ERROR; t <= DIPLO; t++)
      logp_trans[s][t] = -INFINITY;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"\n### [Intvl %d] @ (%d, %d): %d <- %d\n",i,ri.i,ri.j,ri.ci,ri.cj);
#endif

  // Compute and normalize transition probabilities so that \sum_{t'} p(t' <- s) = 1
  psum = 0.;
  for (int s = ERROR; s <= DIPLO; s++)
    { if (dp[i+1][s] == -INFINITY)
        continue;
      next_st = st[i+1][s];
      for (int t = ERROR; t <= DIPLO; t++)
        { if (t == ERROR)
            { logp = logp_e_bw(i,rintvl,plen,perror,cov);
            }
          else if (t == REPEAT)
            { // prev H < this R
              icms = minus_sigma(ic,N_SIGMA_FLUC);
              nc = next_st[HAPLO][0];
              nb = next_st[HAPLO][1];
              logp_fluc = logp_fluctuation(ri.j-1,nb,icms,nc,cov[HAPLO]);
              hr_violated = !(nb == POS_INIT) && !(nc < icms) && !(logp_fluc > THRES_LOGP_FLUC);

              // prev D < this R
              nc = next_st[DIPLO][0];
              nb = next_st[DIPLO][1];
              logp_fluc = logp_fluctuation(ri.j-1,nb,icms,nc,cov[DIPLO]);
              dr_violated = !(nb == POS_INIT) && !(nc < icms) && !(logp_fluc > THRES_LOGP_FLUC);

              if (hr_violated || dr_violated)
                logp = -INFINITY;
              else
                logp = logp_r_bw(i,rintvl,cov,next_st[REPEAT]);
            }
          else if (t == HAPLO)
            { // this H < prev D
              icps = plus_sigma(ic,N_SIGMA_FLUC);
              nc = next_st[DIPLO][0];
              nb = next_st[DIPLO][1];
              logp_fluc = logp_fluctuation(ri.j-1,nb,icps,nc,cov[DIPLO]);
              hd_violated = !(nb == POS_INIT) && !(icps < nc) && !(logp_fluc > THRES_LOGP_FLUC);

              if (hd_violated)
                logp = -INFINITY;
              else
                logp = logp_h_bw(i,rintvl,cerror,cov,next_st[HAPLO]);
            }
          else   // DIPLO
            { // prev H < this D
              icms = minus_sigma(ic,N_SIGMA_FLUC);
              nc = next_st[HAPLO][0];
              nb = next_st[HAPLO][1];
              logp_fluc = logp_fluctuation(ri.j-1,nb,icms,nc,cov[HAPLO]);
              hd_violated = !(nb == POS_INIT) && !(nc < icms) && !(logp_fluc > THRES_LOGP_FLUC);

              // this H < prev D
              icps = plus_sigma(ic,N_SIGMA_FLUC);
              nc = next_st[REPEAT][0];
              nb = next_st[REPEAT][1];
              logp_fluc = logp_fluctuation(ri.j-1,nb,icps,nc,cov[REPEAT]);
              dr_violated = !(nb == POS_INIT) && !(icps < nc) && !(logp_fluc > THRES_LOGP_FLUC);

              if (hd_violated || dr_violated)
                logp = -INFINITY;
              else
                logp = logp_d_bw(i,rintvl,cerror,cov,next_st[DIPLO]);
            }

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
          fprintf(stderr,"%c <- %c: logp = %.lf\n",stoc[t],stoc[s],logp);
#endif

          logp_trans[s][t] = exp(logp);
          psum += logp_trans[s][t];
        }
    }
  for (int s = ERROR; s <= DIPLO; s++)
    { if (dp[i+1][s] == -INFINITY)
        continue;
      for (int t = ERROR; t <= DIPLO; t++)
        logp_trans[s][t] = (psum > 0.) ? log(logp_trans[s][t]/psum) : -INFINITY;
    }

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"Normalized logp_trans:\n");
//   for (int s = ERROR; s <= DIPLO; s++)
//     for (int t = ERROR; t <= DIPLO; t++)
//       fprintf(stderr,"  %c <- %c: logp = %.lf\n",stoc[t],stoc[s],logp_trans[s][t]);
// #endif

  bool only_r = true;
  for (int t = ERROR; t <= DIPLO; t++)
    { max_logp = -INFINITY;
      max_s = -1;
      for (int s = ERROR; s <= DIPLO; s++)
        { if (dp[i+1][s] == -INFINITY)
            continue;
          logp = dp[i+1][s]+logp_trans[s][t];
          if (max_logp < logp)
            { max_logp = logp;
              max_s = s;
              if (t != REPEAT && logp_trans[s][t] > THRES_LOGP_REP)
                only_r = false;
            }
        }
    }

  if (only_r)
    { 
#ifdef DEBUG_DP
      fprintf(stderr,"Absolutely REPEAT @ %d\n",i);
#endif
      rep[i] = true;
      for (int s = ERROR; s <= DIPLO; s++)
        { dp[i][s] = dp[i+1][s];
          if (dp[i+1][s] == -INFINITY)
            continue;
          for (int ii = i+1; ii < M; ii++)
            bt[i][s][ii] = bt[i+1][s][ii];
          bt[i][s][i] = s;
          for (int t = ERROR; t <= DIPLO; t++)
            for (int n = 0; n < 2; n++)
              st[i][s][t][n] = st[i+1][s][t][n];
        }
      rintvl[i].i = rintvl[i+1].i;
      rintvl[i].j = rintvl[i+1].j;
      rintvl[i].ci = rintvl[i+1].ci;
      rintvl[i].cj = rintvl[i+1].cj;
      rintvl[i].asgn = rintvl[i+1].asgn;
      return;
    }

  int curr_h, curr_d, curr_r;
  double dh_ratio;

  for (int t = ERROR; t <= DIPLO; t++)
    { max_logp = -INFINITY;
      max_s = -1;
      for (int s = ERROR; s <= DIPLO; s++)
        { if (dp[i+1][s] == -INFINITY)
          continue;
        logp = dp[i+1][s]+logp_trans[s][t];
        if (max_logp < logp)
          { max_logp = logp;
            max_s = s;
          }
        }
      dp[i][t] = max_logp;
      if (max_s != -1)
        { for (int ii = i+1; ii < M; ii++)
            bt[i][t][ii] = bt[i+1][max_s][ii];
          bt[i][t][i] = t;

          if (t == ERROR)
            { for (int s = REPEAT; s <= DIPLO; s++)
                for (int n = 0; n < 2; n++)
                  st[i][t][s][n] = st[i+1][max_s][s][n];
            }
          else if (t == REPEAT)
            { for (int s = HAPLO; s <= DIPLO; s++)
                for (int n = 0; n < 2; n++)
                  st[i][t][s][n] = st[i+1][max_s][s][n];
              if (st[i+1][max_s][REPEAT][0] < MIN(cov[REPEAT],ri.ci))
                for (int n = 0; n < 2; n++)
                  st[i][t][REPEAT][n] = st[i+1][max_s][REPEAT][n];
              else
                { st[i][t][REPEAT][0] = MIN(cov[REPEAT],ri.ci);
                  st[i][t][REPEAT][1] = ri.i;
                }
            }
          else if (t == HAPLO)
            { curr_h = ri.ci;
              dh_ratio = calc_h_to_d_bw(i,bt[i][t],rintvl,M);
              if (dh_ratio == -INFINITY)
                { bool has_d = false;
                  for (int ii = i+1; ii < M; ii++)
                    if (bt[i][t][ii] == DIPLO)
                      has_d = true;
                  if (has_d)
                    curr_d = st[i+1][max_s][DIPLO][0];   // TODO: update pos as well?
                  else
                    curr_d = MAX(cov[DIPLO],curr_h*2);
                }
              else
                curr_d = (int)(dh_ratio*curr_h);
              curr_r = (int)(DR_RATIO_R*curr_d);
              st[i][t][HAPLO][0] = curr_h;
              st[i][t][HAPLO][1] = ri.i;
              st[i][t][DIPLO][0] = curr_d;
              st[i][t][DIPLO][1] = ri.i+POS_OFFSET;
              st[i][t][REPEAT][0] = curr_r;
              st[i][t][REPEAT][1] = ri.i+POS_OFFSET;
            }
          else   // DIPLO
            { curr_d = ri.ci;
              dh_ratio = calc_d_to_h_bw(i,bt[i][t],rintvl,M);;
              if (dh_ratio == -INFINITY)
                { bool has_h = false;
                  for (int ii = i+1; ii < M; ii++)
                    if (bt[i][t][ii] == HAPLO)
                      has_h = true;
                  if (has_h)
                    curr_h = st[i+1][max_s][HAPLO][0];   // TODO: update pos as well?
                  else
                    curr_h = curr_d/2;
                }
              else
                curr_h = (int)((double)curr_d/dh_ratio);
              curr_r = (int)(DR_RATIO_R*curr_d);
              st[i][t][HAPLO][0] = curr_h;
              st[i][t][HAPLO][1] = ri.i+POS_OFFSET;
              st[i][t][DIPLO][0] = curr_d;
              st[i][t][DIPLO][1] = ri.i;
              st[i][t][REPEAT][0] = curr_r;
              st[i][t][REPEAT][1] = ri.i+POS_OFFSET;
            }
        }

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      for (int s = ERROR; s <= DIPLO; s++)
        { logp = dp[i+1][s]+logp_trans[s][t];
          fprintf(stderr,"%c (<- %c): logp = %.lf\n",stoc[t],stoc[s],logp);
        }
#endif
    }

  return;
}

static inline void backward_print(int i, double **dp, int ****st, char ***bt, Rel_Intvl *rintvl, int M)
{ Rel_Intvl ri = rintvl[i];
  double max_logp = -INFINITY;
  int max_s = -1;

  for (int s = ERROR; s <= DIPLO; s++)
    { if (max_logp < dp[i][s])
        { max_logp = dp[i][s];
          max_s = s;
        }
    }

  fprintf(stderr,"\n@@@ Intvl[%d] (%d, %d): %d <- %d\n",i,ri.i,ri.j,ri.ci,ri.cj);
  for (int s = ERROR; s <= DIPLO; s++)
    { if (i == M-1)
        { fprintf(stderr,"  %c: logp=%5.lf, H/D/R=(%2d,%2d,%3d) %s\n",stoc[s],dp[i][s],st[i][s][HAPLO][0],st[i][s][DIPLO][0],st[i][s][REPEAT][0],(s == max_s) ? "***" : "");
        }
      else
        { int next_s = bt[i][s][i+1];
          if (next_s == ' ')   // impossible state
            continue;
          int **next_st = st[i+1][next_s];
          fprintf(stderr,"  %c (<-%c; %2d,%2d,%3d): logp=%5.lf, H/D/R=(%2d,%2d,%3d) ",stoc[s],stoc[next_s],next_st[HAPLO][0],next_st[DIPLO][0],next_st[REPEAT][0],dp[i][s],st[i][s][HAPLO][0],st[i][s][DIPLO][0],st[i][s][REPEAT][0]); 
          for (int ii = i; ii < M; ii++)
            fprintf(stderr,"%c",stoc[(int)bt[i][s][ii]]);
          fprintf(stderr," %s\n",(s == max_s) ? "***" : "");
        }
    }
  fprintf(stderr,"\n");
  
  return;
}

static void backward_dp(double **dp, int ****st, char ***bt, bool *rep, char *asgn, Rel_Intvl *_rintvl, int M, int plen, P_Error *perror, P_Error *cerror, int cov[])
{ // Copy rintvl
  Rel_Intvl *rintvl = Malloc(sizeof(Rel_Intvl)*M,"Copy rintvl");
  for (int i = 0; i < M; i++)
    { rintvl[i].i = _rintvl[i].i;
      rintvl[i].j = _rintvl[i].j;
      rintvl[i].ci = _rintvl[i].ci;
      rintvl[i].cj = _rintvl[i].cj;
      rintvl[i].asgn = _rintvl[i].asgn;
    }
  
  // Initialization
  for (int i = 0; i < M; i++)
    { rep[i] = false;
      for (int s = ERROR; s <= DIPLO; s++)
        { dp[i][s] = -INFINITY;
          for (int t = REPEAT; t <= DIPLO; t++)
            st[i][s][t][0] = -1;
          for (int j = 0; j < M; j++)
            bt[i][s][j] = ' ';
        }
    }

  const int POS_INIT = plen+POS_OFFSET;
  int i;

  // First interval
  i = M-1;
  Rel_Intvl ri = rintvl[i];
  for (int s = ERROR; s <= DIPLO; s++)
    { for (int t = REPEAT; t <= DIPLO; t++)
        { st[i][s][t][0] = cov[t];
          st[i][s][t][1] = POS_INIT;
        }
      bt[i][s][i] = s;
    }
  
  dp[i][ERROR] = logp_e_bw(i,rintvl,plen,perror,cov);
  dp[i][REPEAT] = logp_r_bw(i,rintvl,cov,st[i][REPEAT][REPEAT]);
  
  dp[i][HAPLO] = logp_poisson(ri.cj, cov[HAPLO]);
  st[i][HAPLO][HAPLO][0] = ri.ci;
  st[i][HAPLO][HAPLO][1] = ri.i;
  st[i][HAPLO][DIPLO][0] = ri.ci*2;
  st[i][HAPLO][DIPLO][1] = ri.i+POS_OFFSET;
  st[i][HAPLO][REPEAT][0] = (int)(DR_RATIO_R*2*ri.ci);
  st[i][HAPLO][REPEAT][1] = ri.i+POS_OFFSET;
  
  dp[i][DIPLO] = logp_poisson(ri.cj, cov[DIPLO]);
  st[i][DIPLO][HAPLO][0] = ri.ci;
  st[i][DIPLO][HAPLO][1] = ri.i+POS_OFFSET;
  st[i][DIPLO][DIPLO][0] = ri.ci;
  st[i][DIPLO][DIPLO][1] = ri.i;
  st[i][DIPLO][REPEAT][0] = (int)(DR_RATIO_R*ri.ci);
  st[i][DIPLO][REPEAT][1] = ri.i+POS_OFFSET;

  double psum = 0.;
  for (int s = ERROR; s <= DIPLO; s++)
    psum += exp(dp[i][s]);
  for (int s = ERROR; s <= DIPLO; s++)
    dp[i][s] = log(exp(dp[i][s])/psum);   // TODO: reduce log/exp

#ifdef DEBUG_DP
  backward_print(i,dp,st,bt,rintvl,M);
#endif

  // The others
  for (i = M-2; i >= 0; i--)
    { backward_update(i,dp,st,bt,rep,POS_INIT,rintvl,M,plen,perror,cerror,cov);

#ifdef DEBUG_DP
      backward_print(i,dp,st,bt,rintvl,M);
#endif
    }

  // Traceback
  double max_logp;
  int max_s;

  max_logp = -INFINITY;
  max_s = -1;
  for (int s = ERROR; s <= DIPLO; s++)
    { if (max_logp < dp[0][s])
        { max_logp = dp[0][s];
          max_s = s;
        }
    }
  for (i = 0; i < M; i++)
    asgn[i] = (rep[i]) ? REPEAT : bt[0][max_s][i];

  free(rintvl);
  return;
}

/*******************************************************************************************
 *
 *  Main routine
 *
 ********************************************************************************************/

void classify_reliable(Rel_Intvl *rintvl, int M, Intvl *intvl, int N, int plen, P_Error *perror, P_Error *cerror, int hcov, int dcov)
{ int cov[N_STATE] = {1, dcov+6*(int)sqrt(dcov), hcov, dcov};
  
  DR_RATIO_R = 1.+N_SIGMA_R*(1./sqrt(cov[DIPLO]));

  if (M == 0)   // TODO: use global cov?
    return;

  double **dp_f, **dp_b;
  int ****st_f, ****st_b;
  char ***bt_f, ***bt_b;
  bool *rep_f, *rep_b;
  char *asgn_f, *asgn_b;

  dp_f = Malloc(sizeof(double*)*M,"");
  for (int i = 0; i < M; i++)
    dp_f[i] = Malloc(sizeof(double)*N_STATE,"");
  st_f = Malloc(sizeof(int***)*M,"");
  for (int i = 0; i < M; i++)
    { st_f[i] = Malloc(sizeof(int**)*N_STATE,"");
      for (int j = 0; j < N_STATE; j++)
        { st_f[i][j] = Malloc(sizeof(int*)*N_STATE,"");
          for (int k = 0; k < N_STATE; k++)
            st_f[i][j][k] = Malloc(sizeof(int)*2,"");
        }
    }
  bt_f = Malloc(sizeof(char**)*M,"");
  for (int i = 0; i < M; i++)
    { bt_f[i] = Malloc(sizeof(char*)*N_STATE,"");
      for (int j = 0; j < N_STATE; j++)
        bt_f[i][j] = Malloc(sizeof(char)*M,"");
    }
  rep_f = Malloc(sizeof(bool)*M,"");
  asgn_f = Malloc(sizeof(char)*M,"");

  forward_dp(dp_f,st_f,bt_f,rep_f,asgn_f,rintvl,M,plen,perror,cerror,cov);

#ifdef DEBUG_ITER
  fprintf(stderr,"  FWD : ");
  for (int i = 0; i < M; i++)
    fprintf(stderr,"%c",stoc[(int)asgn_f[i]]);
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  dp_b = Malloc(sizeof(double*)*M,"");
  for (int i = 0; i < M; i++)
    dp_b[i] = Malloc(sizeof(double)*N_STATE,"");
  st_b = Malloc(sizeof(int***)*M,"");
  for (int i = 0; i < M; i++)
    { st_b[i] = Malloc(sizeof(int**)*N_STATE,"");
      for (int j = 0; j < N_STATE; j++)
        { st_b[i][j] = Malloc(sizeof(int*)*N_STATE,"");
          for (int k = 0; k < N_STATE; k++)
            st_b[i][j][k] = Malloc(sizeof(int)*2,"");
        }
    }
  bt_b = Malloc(sizeof(char**)*M,"");
  for (int i = 0; i < M; i++)
    { bt_b[i] = Malloc(sizeof(char*)*N_STATE,"");
      for (int j = 0; j < N_STATE; j++)
        bt_b[i][j] = Malloc(sizeof(char)*M,"");
    }
  rep_b = Malloc(sizeof(bool)*M,"");
  asgn_b = Malloc(sizeof(char)*M,"");

  backward_dp(dp_b,st_b,bt_b,rep_b,asgn_b,rintvl,M,plen,perror,cerror,cov);

#ifdef DEBUG_ITER
  fprintf(stderr,"  BWD : ");
  for (int i = 0; i < M; i++)
    fprintf(stderr,"%c",stoc[(int)asgn_b[i]]);
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  // Find best path by combining forward and backward DPs
  double max_logp, logp_f, logp_b, logp_joint;
  int max_s;
  for (int i = 0; i < M; i++)
    { max_logp = -INFINITY;
      max_s = -1;

      for (int s = ERROR; s <= DIPLO; s++)
        { if (rep_f[i])
            logp_f = (s == REPEAT) ? 0. : -INFINITY;
          else
            logp_f = dp_f[i][s];
          if (rep_b[i])
            logp_b = (s == REPEAT) ? 0. : -INFINITY;
          else
            logp_b = dp_b[i][s];
            
          logp_joint = logp_f + logp_b;
          if (max_logp < logp_joint)
            { max_logp = logp_joint;
              max_s = s;
            }
        }
      if (max_s == -1)
        max_s = REPEAT;
      rintvl[i].asgn = max_s;
    }

#ifdef DEBUG_ITER
  fprintf(stderr,"  REL : ");
  for (int i = 0; i < M; i++)
    fprintf(stderr,"%c",stoc[(int)rintvl[i].asgn]);
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

  for (int i = 0; i < M; i++)
    free(dp_f[i]);
  free(dp_f);
  for (int i = 0; i < M; i++)
    { for (int j = 0; j < N_STATE; j++)
        { for (int k = 0; k < N_STATE; k++)
            free(st_f[i][j][k]);
          free(st_f[i][j]);
        }
      free(st_f[i]);
    }
  free(st_f);
  for (int i = 0; i < M; i++)
    { for (int j = 0; j < N_STATE; j++)
        free(bt_f[i][j]);
      free(bt_f[i]);
    }
  free(bt_f);
  free(rep_f);
  free(asgn_f);

  for (int i = 0; i < M; i++)
    free(dp_b[i]);
  free(dp_b);
  for (int i = 0; i < M; i++)
    { for (int j = 0; j < N_STATE; j++)
        { for (int k = 0; k < N_STATE; k++)
            free(st_b[i][j][k]);
          free(st_b[i][j]);
        }
      free(st_b[i]);
    }
  free(st_b);
  for (int i = 0; i < M; i++)
    { for (int j = 0; j < N_STATE; j++)
        free(bt_b[i][j]);
      free(bt_b[i]);
    }
  free(bt_b);
  free(rep_b);
  free(asgn_b);

  return;
}
