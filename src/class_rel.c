#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "ClassPro.h"

double DR_RATIO;

Rel_Arg *alloc_rel_arg(int rlen_max)
{ Rel_Arg *arg = Malloc(sizeof(Rel_Arg),"Allocating rel arg");
  arg->COV     = Malloc(N_STATE*sizeof(cnt_t),"COV array");
  arg->dp      = Malloc((rlen_max*N_STATE)*sizeof(double),"DP array");
  arg->st      = Malloc((rlen_max*N_STATE)*sizeof(CovsT),"Cov state array");
  arg->bt      = Malloc((rlen_max*N_STATE)*sizeof(char*),"State path array");
  for (int i = 0; i < rlen_max*N_STATE; i++)
    arg->bt[i] = Malloc(rlen_max*sizeof(char),"State path");
  arg->dh_ratio = Malloc((rlen_max*N_STATE)*sizeof(double),"D/H ratio array");
  arg->rpos     = Malloc(rlen_max*sizeof(bool),"Absolute R array");
  arg->intvl    = Malloc(rlen_max*sizeof(Intvl),"Copied rel intvl array");
  return arg;
}

void free_rel_arg(Rel_Arg *arg, int rlen_max)
{ free(arg->COV);
  free(arg->dp);
  free(arg->st);
  for (int i = 0; i < rlen_max*N_STATE; i++)
    free(arg->bt[i]);
  free(arg->bt);
  free(arg->dh_ratio);
  free(arg->rpos);
  free(arg->intvl);
  free(arg);
  return;
}

static inline int _pred(int x, bool FORWARD)
{ return (FORWARD) ? x-1 : x+1; }

static inline int _succ(int x, bool FORWARD)
{ return (FORWARD) ? x+1 : x-1; }

static inline pos_t _offset(pos_t x, bool FORWARD)
{ return (FORWARD) ? x - OFFSET : x + OFFSET; }

static inline pos_t _beg_pos(Intvl I, bool FORWARD)
{ return (FORWARD) ? I.b : I.e-1; }

static inline cnt_t _beg_cnt(Intvl I, bool FORWARD)
{ return (FORWARD) ? I.ccb : I.cce; }

static inline pos_t _end_pos(Intvl I, bool FORWARD)
{ return (FORWARD) ? I.e-1 : I.b; }

static inline cnt_t _end_cnt(Intvl I, bool FORWARD)
{ return (FORWARD) ? I.cce : I.ccb; }

static inline enum State _find_max_dp(double *dp, int i)
{ double max_logp = -INFINITY;
  int    max_s    = N_STATE;
  for (int s = ERROR; s <= DIPLO; s++)
    { int idx = REL_IDX(i,s);
      if (max_logp < dp[idx])
        { max_logp = dp[idx];
          max_s = s;
        }
    }
  return max_s;   // NOTE: N_STATE means all logp = -inf (impossible state)
}

typedef struct
  { enum State max_x;
    double     max_logp;
  } Max_Cell;

static inline Max_Cell _find_max_dp_tr(double *dp, double logp_tr[N_STATE][N_STATE],
                                       int i, enum State s, enum State t, bool FORWARD)
{ int i_pred = _pred(i, FORWARD);
  double max_logp = -INFINITY;
  int    max_x    = N_STATE;
  for (enum State x = ERROR; x <= DIPLO; x++)
    { enum State _s = (s < N_STATE) ? s : x;
      enum State _t = (t < N_STATE) ? t : x;
      double logp = dp[REL_IDX(i_pred,_s)]+logp_tr[_s][_t];
      if (max_logp < logp)
        { max_logp = logp;
          max_x = x;
        }
    }
  Max_Cell ret = { max_x, max_logp };
  return ret;
}

static inline int find_nn(bool forward, int i, enum State s, char *asgn, int L)
{ int idx = i;
  if (forward)
    while (idx < L && asgn[idx] != (char)s)
      idx++;
  else
    while (idx >= 0 && asgn[idx] != (char)s)
      idx--;
  return idx;
}

static inline bool _is_out(int idx, bool _F, int L)
{ return (_F && idx < 0) || (!_F && idx >= L); }

// L == len(asgn) == len(intvl)
static inline double calc_dh_ratio(enum State init_s, char *asgn, Intvl *intvl, int L, bool FORWARD)
{
// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"[DH_RATIO] L=%d, init_s=%c\n",L,stoc[init_s]);
//   for (int i = 0; i < L; i++)
//     fprintf(stderr,"i = %d: intvl=(%d, %d) asgn=%c\n",i,intvl[i].b,intvl[i].e,stoc[(int)asgn[i]]);
// #endif
          
  int idx[3+1];
  idx[0] = (FORWARD) ? L : -1;
  enum State s = init_s;
  for (int i = 0; i < 3; i++)
    { idx[i+1] = find_nn(!FORWARD,_pred(idx[i], FORWARD),s,asgn,L);
      if (_is_out(idx[i+1],FORWARD,L))
        return -INFINITY;
      s = (s == HAPLO) ? DIPLO : HAPLO;
    }

  Pos_Cnt s1 = { _beg_pos(intvl[idx[1]], FORWARD), _beg_cnt(intvl[idx[1]], FORWARD) };
  Pos_Cnt t  = { _end_pos(intvl[idx[2]], FORWARD), _end_cnt(intvl[idx[2]], FORWARD) };
  Pos_Cnt s2 = { _end_pos(intvl[idx[3]], FORWARD), _end_cnt(intvl[idx[3]], FORWARD) };
  if (!FORWARD)
    { Pos_Cnt tmp = s1;
      s1 = s2;
      s2 = tmp;
    }

// #if defined(DEBUG_REL) && defined(DEBUG_PROB)
//   fprintf(stderr,"L=%d, init_s=%c\n",L,stoc[init_s]);
//   for (int i = 0; i < L; i++)
//     fprintf(stderr,"i = %d: intvl=(%d, %d) asgn=%c\n",i,intvl[i].b,intvl[i].e,stoc[(int)asgn[i]]);
//   fprintf(stderr,"  s1 = %d @ %d\n",s1.cnt,s1.pos);
//   fprintf(stderr,"  t  = %d @ %d\n",t.cnt,t.pos);
//   fprintf(stderr,"  s2 = %d @ %d\n",s2.cnt,s2.pos);
// #endif

  double est_s_cnt = linear_interpolation(t.pos,s2,s1);
  double dh_ratio = (init_s == DIPLO) ? est_s_cnt/t.cnt : t.cnt/est_s_cnt;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"D/H ratio=%5.2lf\n",dh_ratio);
#endif
  return dh_ratio;
}

static inline double logp_e(int idx, Intvl *intvl, cnt_t *COV)
{ Intvl I = intvl[idx];
  double logp, logp_po, logp_er;

  logp_er = I.pe;
  logp_po = logp_poisson(I.ccb,COV[ERROR])+logp_poisson(I.cce,COV[ERROR])+E_PO_BASE;
  logp = MAX(logp_po,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [E] logp(PO) = %.1lf, logp(ER) = %.1lf\n",logp_po,logp_er);
#endif
  return logp;
}

static inline double logp_r(int idx, Intvl *intvl, Pos_Cnt st_pred_r, bool FORWARD, cnt_t *COV)
{ Intvl I = intvl[idx];
  // pos_t beg_pos = _beg_pos(I, FORWARD);
  cnt_t beg_cnt = _beg_cnt(I, FORWARD);

  double logp, logp_sf, logp_er;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [R] ");
#endif

  // double _lambda = (double)pc*(ri.i-pe)/READ_LEN;
  // logp_sf = logp_skellam(ri.ci-pc,_lambda);
  logp_sf = -INFINITY;
  logp_er = (beg_cnt < st_pred_r.cnt) ? logp_binom(beg_cnt,st_pred_r.cnt,1-PE_MEAN) : -INFINITY;
  logp = MAX(logp_sf,logp_er);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"logp(SF) = %.1lf, logp(ER) = %.1lf\n",logp_sf,logp_er);
#endif

  if (logp > R_LOGP)
    return logp;
  cnt_t max_cc = MAX(I.ccb,I.cce);
  if (max_cc >= COV[REPEAT])
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,">= global R-cov\n");
#endif
      return R_LOGP;
    }
  if (max_cc >= st_pred_r.cnt)
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,">= est R-cov\n");
#endif
      return R_LOGP;
    }
  return logp;
}

static inline double logp_h(int idx, Intvl *intvl, double *dh_ratio, enum State s, CovsT st_pred, bool FORWARD) 
{ Intvl I = intvl[idx];
  pos_t beg_pos = _beg_pos(I, FORWARD);
  cnt_t beg_cnt = _beg_cnt(I, FORWARD);

  Pos_Cnt st;
  double logp, logp_sf_h, logp_sf_d;

  st = st_pred[HAPLO];
  // fprintf(stderr,"H trans @ H: %d @ %d ~~ %d @ %d (cov = %d)\n",
  //                st.cnt,_pred(st.pos),beg_cnt,beg_pos,st.cnt);
  logp_sf_h = logp_trans(_pred(st.pos, FORWARD),beg_pos,st.cnt,beg_cnt,st.cnt);   // TODO: cache?

  logp_sf_d = 0.;
  double r = dh_ratio[REL_IDX(_pred(idx, FORWARD),s)];
  if (r != -INFINITY)
    { st = st_pred[DIPLO];
      // fprintf(stderr,"D trans @ H: %d @ %d ~~ %d @ %d (r = %lf, cov = %d)\n",
      //                st.cnt,_pred(st.pos),(int)(r*beg_cnt),beg_pos,r,st.cnt);
      logp_sf_h = logp_trans(_pred(st.pos, FORWARD),beg_pos,st.cnt,(int)(r*beg_cnt),st.cnt);
    }
  logp = logp_sf_h+logp_sf_d;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [H] logp = %.1lf + %.1lf = %.1lf\n",logp_sf_h,logp_sf_d,logp);
#endif
  return logp;
}

static inline double logp_d(int idx, Intvl *intvl, double *dh_ratio, enum State s, CovsT st_pred, bool FORWARD) 
{ Intvl I = intvl[idx];
  pos_t beg_pos = _beg_pos(I, FORWARD);
  cnt_t beg_cnt = _beg_cnt(I, FORWARD);

  Pos_Cnt st;
  double logp, logp_sf_h, logp_sf_d;

  logp_sf_d = 0.;
  double r = dh_ratio[REL_IDX(_pred(idx, FORWARD),s)];
  if (r != -INFINITY)
    { st = st_pred[HAPLO];
      // fprintf(stderr,"D trans @ H: %d @ %d ~~ %d @ %d (r = %lf, cov = %d)\n",
      //                st.cnt,_pred(st.pos),(int)((double)beg_cnt/r),beg_pos,r,st.cnt);
      logp_sf_h = logp_trans(_pred(st.pos, FORWARD),beg_pos,st.cnt,(int)((double)beg_cnt/r),st.cnt);
    }

  st = st_pred[DIPLO];
  // fprintf(stderr,"D trans @ D: %d @ %d ~~ %d @ %d (cov = %d)\n",
  //                 st.cnt,_pred(st.pos),beg_cnt,beg_pos,st.cnt);
  logp_sf_h = logp_trans(_pred(st.pos, FORWARD),beg_pos,st.cnt,beg_cnt,st.cnt);

  logp = logp_sf_h+logp_sf_d;

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [D] logp = %.1lf + %.1lf = %.1lf\n",logp_sf_h,logp_sf_d,logp);
#endif
  return logp;
}

static inline double calc_logp(enum State s, enum State t, int idx, Intvl *intvl, CovsT st_pred, double *dh_ratio, bool FORWARD, cnt_t *COV)
{ if (t == ERROR) return logp_e(idx,intvl, COV);
  else if (t == HAPLO) return logp_h(idx,intvl,dh_ratio,s,st_pred, FORWARD);
  else if (t == DIPLO) return logp_d(idx,intvl,dh_ratio,s,st_pred, FORWARD);
  else return logp_r(idx,intvl,st_pred[REPEAT], FORWARD, COV);
}

static void _update(Rel_Arg *arg, int i, int M)
{ bool     FORWARD = arg->FORWARD;
  cnt_t   *COV = arg->COV;
  double  *dp = arg->dp;
  CovsT   *st = arg->st;
  char   **bt = arg->bt;
  double  *dh_ratio = arg->dh_ratio;
  bool    *rpos = arg->rpos;
  Intvl   *intvl = arg->intvl;

  Intvl I = intvl[i];
  pos_t end_pos = _end_pos(I, FORWARD);
  cnt_t end_cnt = _end_cnt(I, FORWARD);
  int i_pred = _pred(i, FORWARD);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"\n### [Intvl %d] @ (%d, %d): %d ~~ %d\n",i,I.b,I.e,I.ccb,I.cce);
  fprintf(stderr,"Raw logp_tr:\n");
#endif

  // Compute and normalize transition probabilities so that \sum_{t'} p(s -> t') = 1
  double logp_tr[N_STATE][N_STATE];   // s @ i_pred -> t @ i
  for (int s = ERROR; s <= DIPLO; s++)
    for (int t = ERROR; t <= DIPLO; t++)
      logp_tr[s][t] = -INFINITY;

  for (int s = ERROR; s <= DIPLO; s++)
    { 
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,"s = %c ->\n",stoc[s]);
#endif
      int idx = REL_IDX(i_pred,s);
      if (dp[idx] == -INFINITY)
        { for (int t = ERROR; t <= DIPLO; t++)
          logp_tr[s][t] = 0.;
          continue;
        }
      for (int t = ERROR; t <= DIPLO; t++)
        { double logp = calc_logp(s,t,i,intvl,st[idx],dh_ratio, FORWARD, COV);
          logp_tr[s][t] = exp(logp);
        }
    }
  double psum = 0.;
  for (int s = ERROR; s <= DIPLO; s++)
    for (int t = ERROR; t <= DIPLO; t++)
      psum += logp_tr[s][t];
#ifdef DEBUG
  if (psum == 0.)
    { fprintf(stderr,"No possible state @ %d\n",i);
      // exit(1);
      for (int s = ERROR; s <= DIPLO; s++)
        logp_tr[s][ERROR] = 1.;
      psum = 4.;
    }
#endif
  for (int s = ERROR; s <= DIPLO; s++)
    for (int t = ERROR; t <= DIPLO; t++)
      logp_tr[s][t] = log(logp_tr[s][t]/psum);

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"Normalized logp_tr:\n");
  for (int t = ERROR; t <= DIPLO; t++)
    { fprintf(stderr,"t = %c:",stoc[t]);
      for (int s = ERROR; s <= DIPLO; s++)
        fprintf(stderr,"   %5.lf (%c)",logp_tr[s][t],stoc[s]);
      fprintf(stderr,"\n");
    }
#endif

  // All paths converge to R at this interval
  bool only_r = true;
  for (int s = ERROR; s <= DIPLO; s++)
    { int maxt = _find_max_dp_tr(dp,logp_tr,i,s,N_STATE,FORWARD).max_x;
      if (maxt != N_STATE && maxt != REPEAT)
        { only_r = false;
          break;
        }
    }
  if (only_r)
    { rpos[i] = true;
      intvl[i] = intvl[i_pred];
      for (int s = ERROR; s <= DIPLO; s++)
        { int idx = REL_IDX(i,s);
          int idx_pred = REL_IDX(i_pred,s);
          dp[idx] = dp[idx_pred];
          if (dp[idx] == -INFINITY)
            continue;
          if (FORWARD)
            for (int ii = 0; ii < i; ii++)
              bt[idx][ii] = bt[idx_pred][ii];
          else
            for (int ii = i+1; ii < M; ii++)
              bt[idx][ii] = bt[idx_pred][ii];
          bt[idx][i] = s;
          for (int t = ERROR; t <= DIPLO; t++)
            st[idx][t] = st[idx_pred][t];
        }
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      fprintf(stderr,"Absolutely REPEAT @ %d\n",i);
#endif
      return;
    }

  // Update
  int maxs_h = _find_max_dp_tr(dp,logp_tr,i,N_STATE,HAPLO,FORWARD).max_x;
  int maxs_d = _find_max_dp_tr(dp,logp_tr,i,N_STATE,DIPLO,FORWARD).max_x;
  if (maxs_h == HAPLO && maxs_d == DIPLO)
    logp_tr[HAPLO][HAPLO] = logp_tr[DIPLO][DIPLO] = MIN(logp_tr[HAPLO][HAPLO],logp_tr[DIPLO][DIPLO]);

  int curr_h, curr_d, curr_r;
  double r;
  for (int t = ERROR; t <= DIPLO; t++)
    { Max_Cell max_cell = _find_max_dp_tr(dp,logp_tr,i,N_STATE,t,FORWARD);
      int max_s = max_cell.max_x;
      double max_logp = max_cell.max_logp;

      int idx = REL_IDX(i,t);
      int idx_pred = REL_IDX(i_pred,max_s);
      dp[idx] = max_logp;
      if (max_s == N_STATE)
        continue;

      if (FORWARD)
        for (int ii = 0; ii < i; ii++)
          bt[idx][ii] = bt[idx_pred][ii];
      else
        for (int ii = i+1; ii < M; ii++)
          bt[idx][ii] = bt[idx_pred][ii];
      bt[idx][i] = t;

      if (t == ERROR)
        { for (int s = REPEAT; s <= DIPLO; s++)
            st[idx][s] = st[idx_pred][s];
        }
      else if (t == REPEAT)
        { for (int s = HAPLO; s <= DIPLO; s++)
            { st[idx][s].pos = _offset(end_pos, FORWARD);
              st[idx][s].cnt = st[idx_pred][s].cnt;
            }
          cnt_t r_cnt = MIN(end_cnt,COV[REPEAT]);
          if (st[idx_pred][REPEAT].cnt < r_cnt)
            st[idx][REPEAT] = st[idx_pred][REPEAT];
          else
            { st[idx][REPEAT].pos = _offset(end_pos, FORWARD);
              st[idx][REPEAT].cnt = r_cnt;
            }
        }
      else if (t == HAPLO)
        { curr_h = end_cnt;
          r = calc_dh_ratio(HAPLO,(FORWARD) ? bt[idx] : bt[idx]+i,
                                  (FORWARD) ? intvl : intvl+i,
                                  (FORWARD) ? i+1 : M-i, FORWARD);
          if (r == -INFINITY)
            { bool has_d = false;
              if (FORWARD)
                { for (int ii = 0; ii < i; ii++)
                    if (bt[idx][ii] == DIPLO)
                      has_d = true;
                }
              else
                { for (int ii = i+1; ii < M; ii++)
                    if (bt[idx][ii] == DIPLO)
                      has_d = true;
                }
              if (has_d)
                curr_d = st[idx_pred][DIPLO].cnt;
              else
                curr_d = curr_h+COV[HAPLO];
            }
          else
            { curr_d = (int)(r*curr_h);
              dh_ratio[idx] = r;
            }
          curr_r = (int)(DR_RATIO*curr_d);
          st[idx][HAPLO].pos = _offset(end_pos, FORWARD);
          st[idx][HAPLO].cnt = curr_h;
          st[idx][DIPLO].pos = _offset(end_pos, FORWARD);
          st[idx][DIPLO].cnt = curr_d;
          st[idx][REPEAT].pos = _offset(end_pos, FORWARD);
          st[idx][REPEAT].cnt = curr_r;
        }
      else   // DIPLO
        { curr_d = end_cnt;
          r = calc_dh_ratio(DIPLO,(FORWARD) ? bt[idx] : bt[idx]+i,
                                  (FORWARD) ? intvl : intvl+i,
                                  (FORWARD) ? i+1 : M-i, FORWARD);
          if (r == -INFINITY)
            { bool has_h = false;
              if (FORWARD)
                { for (int ii = 0; ii < i; ii++)
                    if (bt[idx][ii] == HAPLO)
                      has_h = true;
                }
              else
                { for (int ii = i+1; ii < M; ii++)
                    if (bt[idx][ii] == HAPLO)
                      has_h = true;
                }
              if (has_h)
                curr_h = st[idx_pred][HAPLO].cnt;
              else
                curr_h = MAX(curr_d/2,curr_d-COV[HAPLO]);
            }
          else
            { curr_h = (int)((double)curr_d/r);
              dh_ratio[idx] = r;
            }
          curr_r = (int)(DR_RATIO*curr_d);
          st[idx][HAPLO].pos = _offset(end_pos, FORWARD);
          st[idx][HAPLO].cnt = curr_h;
          st[idx][DIPLO].pos = _offset(end_pos, FORWARD);
          st[idx][DIPLO].cnt = curr_d;
          st[idx][REPEAT].pos = _offset(end_pos, FORWARD);
          st[idx][REPEAT].cnt = curr_r;
        }

      // H < D < R
      if (!( (st[idx][HAPLO].cnt < st[idx][DIPLO].cnt)
             && (st[idx][DIPLO].cnt < st[idx][REPEAT].cnt) ))
        dp[idx] = -INFINITY;
    }

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  for (int s = ERROR; s <= DIPLO; s++)
    { fprintf(stderr,"logp(%c) = %5.lf",stoc[s],dp[REL_IDX(i,s)]);
      for (int t = ERROR; t <= DIPLO; t++)
        fprintf(stderr,"   %c-cov=%3d",stoc[t],st[REL_IDX(i,s)][t].cnt);
      fprintf(stderr,"\n");
    }
  int maxs = _find_max_dp(dp,i);
  fprintf(stderr,"Best state = %c\n",stoc[maxs]);
#endif

  return;
}

char *_classify_rel(Rel_Arg *arg, Intvl *rintvl, int M, int plen)
{ bool     FORWARD = arg->FORWARD;
  cnt_t   *COV = arg->COV;
  double  *dp = arg->dp;
  CovsT   *st = arg->st;
  char   **bt = arg->bt;
  double  *dh_ratio = arg->dh_ratio;
  bool    *rpos = arg->rpos;
  Intvl   *intvl = arg->intvl;

  int idx;
  for (int i = 0; i < M; i++)
    { for (int s = ERROR; s <= DIPLO; s++)
        { idx = REL_IDX(i,s);
          dp[idx] = -INFINITY;
          dh_ratio[idx] = -INFINITY;
          // for (int t = REPEAT; t <= DIPLO; t++)
          //   { st[idx][t].cnt = 0;
          //   }
          // for (int j = 0; j < M; j++)
          //   bt[idx][j] = '';
        }
      rpos[i] = false;
      intvl[i] = rintvl[i];
    }

  const int POS_INIT = _offset((FORWARD) ? 0 : plen, FORWARD);
  int i = (FORWARD) ? 0 : M-1;

  // 1. Init
  Intvl I = intvl[i];
#if defined(DEBUG_REL) && defined(DEBUG_PROB)
  fprintf(stderr,"\n### [Intvl %d] @ (%d, %d): %d ~~ %d\n",i,I.b,I.e,I.ccb,I.cce);
#endif
  for (int s = ERROR; s <= DIPLO; s++)
    { idx = REL_IDX(i,s);
      for (int t = REPEAT; t <= DIPLO; t++)
        { st[idx][t].pos = POS_INIT;
          st[idx][t].cnt = COV[t];
        }
      bt[idx][i] = s;
    }

  idx = REL_IDX(i,ERROR);
  dp[idx] = logp_e(i,intvl, COV);
  
  idx = REL_IDX(i,REPEAT);
  dp[idx] = logp_r(i,intvl,st[idx][REPEAT], FORWARD, COV);
  st[idx][REPEAT].pos = _end_pos(I, FORWARD);
  st[idx][REPEAT].cnt = MIN(_end_cnt(I, FORWARD),COV[REPEAT]);

  idx = REL_IDX(i,HAPLO);
  dp[idx] = logp_poisson(_beg_cnt(I, FORWARD),COV[HAPLO]);
  st[idx][HAPLO].pos = _end_pos(I, FORWARD);
  st[idx][HAPLO].cnt = _end_cnt(I, FORWARD);
  st[idx][DIPLO].pos = _offset(_end_pos(I, FORWARD), FORWARD);
  st[idx][DIPLO].cnt = _end_cnt(I, FORWARD)+COV[HAPLO];

  idx = REL_IDX(i,DIPLO);
  dp[idx] = logp_poisson(_beg_cnt(I, FORWARD), COV[DIPLO]);
  st[idx][HAPLO].pos = _offset(_end_pos(I, FORWARD), FORWARD);
  st[idx][HAPLO].cnt = MAX(_end_cnt(I, FORWARD)/2,(int)_end_cnt(I, FORWARD)-COV[HAPLO]);
  st[idx][DIPLO].pos = _end_pos(I, FORWARD);
  st[idx][DIPLO].cnt = _end_cnt(I, FORWARD);

  double psum = 0.;
  for (int s = ERROR; s <= DIPLO; s++)
    psum += exp(dp[REL_IDX(i,s)]);
  for (int s = ERROR; s <= DIPLO; s++)
    dp[REL_IDX(i,s)] = log(exp(dp[REL_IDX(i,s)])/psum);   // TODO: reduce log/exp if possible

#if defined(DEBUG_REL) && defined(DEBUG_PROB)
      for (int s = ERROR; s <= DIPLO; s++)
        { fprintf(stderr,"logp(%c) = %5.lf",stoc[s],dp[REL_IDX(i,s)]);
          for (int t = ERROR; t <= DIPLO; t++)
            fprintf(stderr,"   %c-cov=%3d",stoc[t],st[REL_IDX(i,s)][t].cnt);
          fprintf(stderr,"\n");
        }
      int maxs = _find_max_dp(dp,i);
      fprintf(stderr,"Best state = %c\n",stoc[maxs]);
#endif

  // 2. Update
  while (true)
    { i = _succ(i, FORWARD);
      if ((FORWARD && i >= M) || (!FORWARD && i < 0))
        break;
      _update(arg,i,M);
    }

  // 3. Traceback
  i = (FORWARD) ? M-1 : 0;
  enum State max_s = _find_max_dp(dp,i);
  idx = REL_IDX(i,max_s);
  for (int j = 0; j < M; j++)
    if (rpos[j])
      bt[idx][j] = REPEAT;

  return bt[idx];
}

typedef struct
  { char  *asgn;
    int    d_diff;
    int    h_diff;
    double hdrr;
  } Iter_Rel;

Iter_Rel classify_rel_fw(Rel_Arg *arg, Intvl *rintvl, int M, int plen)
{ cnt_t *COV = arg->COV;
  arg->FORWARD = true;
  for (int s = ERROR; s <= DIPLO; s++)
    COV[s] = GLOBAL_COV[s];
  char *asgn = _classify_rel(arg,rintvl,M,plen);
  bool no_h = true;
  for (int i = 0; i < M; i++)
    if (asgn[i] == HAPLO)
      no_h = false;
  if (no_h)
    { int l, lsum = 0, csum = 0;
      int first_d_idx = -1;
      for (int i = 0; i < M; i++)
        { if (asgn[i] == DIPLO)
            { l = rintvl[i].e-rintvl[i].b;
              lsum += l;
              csum += (rintvl[i].ccb+rintvl[i].cce)*l/2;
              if (first_d_idx == -1)
                first_d_idx = i;
            }
        }
      if (first_d_idx >= 0)
        { double mean_dcov = (double)csum/lsum;
          if (mean_dcov < GLOBAL_COV[DIPLO])
            { COV[HAPLO] = rintvl[first_d_idx].ccb;
              COV[DIPLO] = COV[HAPLO]+GLOBAL_COV[HAPLO];
              asgn = _classify_rel(arg,rintvl,M,plen);
              no_h = true;
              for (int i = 0; i < M; i++)
                if (asgn[i] == HAPLO)
                  no_h = false;
              if (no_h)
                { lsum = 0; csum = 0;
                  for (int i = 0; i < M; i++)
                    { if (asgn[i] == DIPLO)
                        { l = rintvl[i].e-rintvl[i].b;
                          lsum += l;
                          csum += (rintvl[i].ccb+rintvl[i].cce)*l/2;
                        }
                    }
                  mean_dcov = (double)csum/lsum;
                  if (fabs(mean_dcov-GLOBAL_COV[HAPLO]) <= fabs(mean_dcov-GLOBAL_COV[DIPLO]))
                    for (int i = 0; i < M; i++)
                      if (asgn[i] == DIPLO)
                        asgn[i] = HAPLO;
                }
            }
        }
    }
  // TODO: make a function to get 1) count # of H/D intervals and 2) mean cov of H/D intervals
  { bool all_h = true;
    for (int i = 0; i < M; i++)
      if (asgn[i] != HAPLO)
        all_h = false;
    if (all_h)
      { int l, lsum = 0, csum = 0;
        for (int i = 0; i < M; i++)
          { l = rintvl[i].e-rintvl[i].b;
            lsum += l;
            csum += (rintvl[i].ccb+rintvl[i].cce)*l/2;
          }
        double mean_hcov = (double)csum/lsum;
        if (fabs(mean_hcov-GLOBAL_COV[HAPLO]) >= fabs(mean_hcov-GLOBAL_COV[DIPLO]))
          for (int i = 0; i < M; i++)
            asgn[i] = DIPLO;
      }
  }

  { int n = 0;
    for (int i = 0; i < M; i++)
      if (asgn[i] == HAPLO) n++;
    if (n >= M * 0.7)
      { int l, lsum = 0, csum = 0;
        for (int i = 0; i < M; i++)
          { if (asgn[i] == HAPLO)
            { l = rintvl[i].e-rintvl[i].b;
              lsum += l;
              csum += (rintvl[i].ccb+rintvl[i].cce)*l/2;
            }
          }
        double mean_hcov = (double)csum/lsum;
        if (fabs(mean_hcov-GLOBAL_COV[HAPLO]) >= fabs(mean_hcov-GLOBAL_COV[DIPLO]))
          for (int i = 0; i < M; i++)
            { if (asgn[i] == HAPLO)
                asgn[i] = DIPLO;
              else if (asgn[i] == DIPLO)
                asgn[i] = REPEAT;
            }
      }
  }

  int first_d_idx = -1, last_d_idx = -1;
  int first_h_idx = -1, last_h_idx = -1;
  for (int i = 0; i < M; i++)
    { if (asgn[i] == DIPLO)
        { if (first_d_idx == -1)
            first_d_idx = i;
          last_d_idx = i;
        }
      else if (asgn[i] == HAPLO)
        { if (first_h_idx == -1)
            first_h_idx = i;
          last_h_idx = i;
        }
    }
  int d_diff = (first_d_idx >= 0) ? abs(rintvl[first_d_idx].ccb-rintvl[last_d_idx].cce) : 0;
  int h_diff = (first_h_idx >= 0) ? abs(rintvl[first_h_idx].ccb-rintvl[last_h_idx].cce) : 0;
  double hdrr = (first_d_idx >= 0 && first_h_idx >= 0) ? ((double)rintvl[first_d_idx].ccb/rintvl[first_h_idx].ccb)/((double)rintvl[last_d_idx].cce/rintvl[last_h_idx].cce) : 1.;

  Iter_Rel ret = { asgn, d_diff, h_diff, hdrr };
  return ret;
}

Iter_Rel classify_rel_bw(Rel_Arg *arg, Intvl *rintvl, int M, int plen)
{ cnt_t *COV = arg->COV;
  arg->FORWARD = false;
  for (int s = ERROR; s <= DIPLO; s++)
    COV[s] = GLOBAL_COV[s];
  char *asgn = _classify_rel(arg,rintvl,M,plen);
  bool no_h = true;
  for (int i = 0; i < M; i++)
    if (asgn[i] == HAPLO)
      no_h = false;
  if (no_h)
    { int l, lsum = 0, csum = 0;
      int last_d_idx = -1;
      for (int i = 0; i < M; i++)
        { if (asgn[i] == DIPLO)
            { l = rintvl[i].e-rintvl[i].b;
              lsum += l;
              csum += (rintvl[i].ccb+rintvl[i].cce)*l/2;
              last_d_idx = i;
            }
        }
      if (last_d_idx >= 0)
        { double mean_dcov = (double)csum/lsum;
          if (mean_dcov < GLOBAL_COV[DIPLO])
            { COV[HAPLO] = rintvl[last_d_idx].cce;
              COV[DIPLO] = COV[HAPLO]+GLOBAL_COV[HAPLO];
              asgn = _classify_rel(arg,rintvl,M,plen);
              no_h = true;
              for (int i = 0; i < M; i++)
                if (asgn[i] == HAPLO)
                  no_h = false;
              if (no_h)
                { lsum = 0; csum = 0;
                  for (int i = 0; i < M; i++)
                    { if (asgn[i] == DIPLO)
                        { l = rintvl[i].e-rintvl[i].b;
                          lsum += l;
                          csum += (rintvl[i].ccb+rintvl[i].cce)*l/2;
                        }
                    }
                  mean_dcov = (double)csum/lsum;
                  if (fabs(mean_dcov-GLOBAL_COV[HAPLO]) <= fabs(mean_dcov-GLOBAL_COV[DIPLO]))
                    for (int i = 0; i < M; i++)
                      if (asgn[i] == DIPLO)
                        asgn[i] = HAPLO;
                }
            }
        }
    }
  bool all_h = true;
  for (int i = 0; i < M; i++)
    if (asgn[i] != HAPLO)
      all_h = false;
  if (all_h)
    { int l, lsum = 0, csum = 0;
      for (int i = 0; i < M; i++)
        { l = rintvl[i].e-rintvl[i].b;
          lsum += l;
          csum += (rintvl[i].ccb+rintvl[i].cce)*l/2;
        }
      double mean_hcov = (double)csum/lsum;
      if (fabs(mean_hcov-GLOBAL_COV[HAPLO]) >= fabs(mean_hcov-GLOBAL_COV[DIPLO]))
        for (int i = 0; i < M; i++)
          asgn[i] = DIPLO;
    }

  int n = 0;
  for (int i = 0; i < M; i++)
    if (asgn[i] == HAPLO) n++;
  if (n >= M * 0.7)
    { int l, lsum = 0, csum = 0;
      for (int i = 0; i < M; i++)
        { if (asgn[i] == HAPLO)
          { l = rintvl[i].e-rintvl[i].b;
            lsum += l;
            csum += (rintvl[i].ccb+rintvl[i].cce)*l/2;
          }
        }
      double mean_hcov = (double)csum/lsum;
      if (fabs(mean_hcov-GLOBAL_COV[HAPLO]) >= fabs(mean_hcov-GLOBAL_COV[DIPLO]))
        for (int i = 0; i < M; i++)
          { if (asgn[i] == HAPLO)
              asgn[i] = DIPLO;
            else if (asgn[i] == DIPLO)
              asgn[i] = REPEAT;
          }
    }

  int first_d_idx = -1, last_d_idx = -1;
  int first_h_idx = -1, last_h_idx = -1;
  for (int i = 0; i < M; i++)
    { if (asgn[i] == DIPLO)
        { if (first_d_idx == -1)
            first_d_idx = i;
          last_d_idx = i;
        }
      else if (asgn[i] == HAPLO)
        { if (first_h_idx == -1)
            first_h_idx = i;
          last_h_idx = i;
        }
    }
  int d_diff = (first_d_idx >= 0) ? abs(rintvl[first_d_idx].ccb-rintvl[last_d_idx].cce) : 0;
  int h_diff = (first_h_idx >= 0) ? abs(rintvl[first_h_idx].ccb-rintvl[last_h_idx].cce) : 0;
  double hdrr = (first_d_idx >= 0 && first_h_idx >= 0) ? ((double)rintvl[first_d_idx].ccb/rintvl[first_h_idx].ccb)/((double)rintvl[last_d_idx].cce/rintvl[last_h_idx].cce) : 1;

  Iter_Rel ret = { asgn, d_diff, h_diff, hdrr };
  return ret;
}

static inline bool is_eq_prefix(Intvl *rintvl, int M)
{ if (rintvl[0].asgn != true)
    return false;
  int i = 0;
  while (i < M && rintvl[i].asgn) i++;
  while (i < M)
    { if (rintvl[i].asgn) return false;
      i++;
    }
  return true;
}

static inline bool is_eq_suffix(Intvl *rintvl, int M)
{ if (rintvl[M-1].asgn != true)
    return false;
  int i = M-2;
  while (i >= 0 && rintvl[i].asgn) i--;
  while (i >= 0)
    { if (rintvl[i].asgn) return false;
      i--;
    }
  return true;
}

void classify_rel(Rel_Arg *arg, Intvl *rintvl, int M, Intvl *intvl, int N, int plen)
{ if (M == 0)
    return;

  Iter_Rel cr_f = classify_rel_fw(arg,rintvl,M,plen);
#ifdef DEBUG_REL
#ifdef DEBUG_ITER
  fprintf(stderr,"\n");
#endif
  fprintf(stderr,"  FWD : ");
  for (int i = 0; i < M; i++)
    fprintf(stderr,"%c",stoc[(int)cr_f.asgn[i]]);
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  for (int i = 0; i < M; i++)
    rintvl[i].asgn = cr_f.asgn[i];

  Iter_Rel cr_b = classify_rel_bw(arg,rintvl,M,plen);
#ifdef DEBUG_REL
  fprintf(stderr,"  BWD : ");
  for (int i = 0; i < M; i++)
    fprintf(stderr,"%c",stoc[(int)cr_b.asgn[i]]);
  fprintf(stderr,"\n");
  fflush(stderr);

  fprintf(stderr,"hdrr FWD=%lf, BWD=%lf\n",cr_f.hdrr,cr_b.hdrr);
#endif

  // NOTE: using rintvl as `eqs`
  // for (int i = 0; i < M; i++)
  //   rintvl[i].asgn = (cr_f.asgn[i] == cr_b.asgn[i]) ? true : false;
  bool eq = true;
  for (int i = 0; i < M; i++)
    // if (rintvl[i].asgn == false)
    if (rintvl[i].asgn != cr_b.asgn[i])
      { eq = false;
        break;
      }
  if (eq)
    { //printf("eq\n");
      // for (int i = 0; i < M; i++)
        // rintvl[i].asgn = cr_f.asgn[i];
      //   rintvl[i].asgn = cr_b.asgn[i];
    }
  else
    { //printf("not eq\n");
      if (is_eq_prefix(rintvl,M))
        { //for (int i = 0; i < M; i++)
          //  rintvl[i].asgn = cr_f.asgn[i];
        }
      else if (is_eq_suffix(rintvl,M))
        { for (int i = 0; i < M; i++)
            rintvl[i].asgn = cr_b.asgn[i];
        }
      else
        { if (fabs(cr_f.hdrr-1.) <= fabs(cr_b.hdrr-1.))
          // if (cr_f.d_diff+cr_f.h_diff > cr_b.d_diff+cr_b.h_diff)
            { //for (int i = 0; i < M; i++)
              //  rintvl[i].asgn = cr_f.asgn[i];
            }
          else
            { for (int i = 0; i < M; i++)
                rintvl[i].asgn = cr_b.asgn[i];
            }
        }
    }

#ifdef DEBUG_REL
  fprintf(stderr,"  REL : ");
  for (int i = 0; i < M; i++)
    fprintf(stderr,"%c",stoc[(int)rintvl[i].asgn]);
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  // Copy assignments to `intvl`
  for (int ridx = 0, iidx = 0; ridx < M; ridx++, iidx++)
    { while (iidx < N && !intvl[iidx].is_rel)
        iidx++;
#ifdef DEBUG
      if (iidx >= N || rintvl[ridx].b != intvl[iidx].b || rintvl[ridx].e != intvl[iidx].e)
        { fprintf(stderr,"Inconsistent reliable interval (%d,%d) != (%d,%d)\n",
                  rintvl[ridx].b,rintvl[ridx].e,intvl[iidx].b,intvl[iidx].e);
          exit(1);
        }
#endif
      intvl[iidx].asgn = rintvl[ridx].asgn;
    }

  return;
}
