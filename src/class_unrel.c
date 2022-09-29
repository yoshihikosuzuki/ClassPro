#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "ClassPro.h"

static inline void find_nn_u(int idx, enum State s, Intvl *intvl, int N, int ret[2])
{ int l = idx-1;
  while (l >= 0 && !(intvl[l].asgn == (char)s && intvl[l].is_rel))
    l--;
  if (l < 0)
    l = -1;
  ret[0] = l;
  int r = idx+1;
  while (r < N && !(intvl[r].asgn == (char)s && intvl[r].is_rel))
    r++;
  if (r >= N)
    r = -1;
  ret[1] = r;
  return;
}

static inline cnt_t est_cov(pos_t x, int idx, Intvl *intvl, int N, enum State s, bool from_est)
{ int nn_idx[2];
  find_nn_u(idx,s,intvl,N,nn_idx);
  int l = nn_idx[0], r = nn_idx[1];
  if (l != -1 && r != -1)
    { Intvl L = intvl[l], R = intvl[r];
      Pos_Cnt pc1 = { L.e-1, L.cce }, pc2 = { R.b, R.ccb };
      return (cnt_t)linear_interpolation(x,pc1,pc2);
    }
  else if (l != -1)
    { Intvl L = intvl[l];
      return L.cce;
    }
  else if (r != -1)
    { Intvl R = intvl[r];
      return R.ccb;
    }
  if (from_est)
    return 0;
  cnt_t cov = est_cov(x,idx,intvl,N,(s == HAPLO) ? DIPLO : HAPLO,true);
  if (cov > 0)
    return (s == HAPLO) ? cov/2 : cov*2;
  else
    return GLOBAL_COV[s];
}

static inline double logp_e_u(int idx, Intvl *intvl)
{ Intvl I = intvl[idx];
  double logp, logp_po, logp_er;

  logp_er = I.pe;
  logp_po = logp_poisson(I.cb,GLOBAL_COV[ERROR])+logp_poisson(I.ce,GLOBAL_COV[ERROR])+E_PO_BASE;
  logp = MAX(logp_er,logp_po);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [E] logp(PO)=%lf, logp(ER)=%lf\n",logp_po,logp_er);
#endif
  return logp;
}

static inline double logp_r_u(int idx, Intvl *intvl, int N)
{ Intvl I = intvl[idx];

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
      fprintf(stderr,"  [R] ");
#endif

  if (MAX(I.cb,I.ce) >= GLOBAL_COV[REPEAT])
    { 
#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
      fprintf(stderr,">= Global R-cov\n");
#endif
      return 0.;
    }

  int nn_idx[2];
  find_nn_u(idx,DIPLO,intvl,N,nn_idx);
  int l = nn_idx[0], r = nn_idx[1];
  cnt_t dcov_l, dcov_r;
  if (l == -1 && r == -1)
    dcov_l = dcov_r = GLOBAL_COV[DIPLO];
  else if (l == -1)
    dcov_l = dcov_r = intvl[r].cb;
  else if (r == -1)
    dcov_l = dcov_r = intvl[l].ce;
  else
    { dcov_l = intvl[l].ce;
      dcov_r = intvl[r].cb;
    }
  cnt_t rcov_l = (cnt_t)(DR_RATIO*dcov_l);
  cnt_t rcov_r = (cnt_t)(DR_RATIO*dcov_r);
  if (I.cb >= rcov_l || I.ce >= rcov_r)
    { 
#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
      fprintf(stderr,">= Est R-cov\n");
#endif
      return R_LOGP;
    }
  double logp_l = logp_binom(I.cb,rcov_l,1-PE_MEAN);
  double logp_r = logp_binom(I.ce,rcov_r,1-PE_MEAN);
  double logp = logp_l+logp_r;

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
  fprintf(stderr,"  [R] logp=%lf\n",logp);
#endif
  return logp;
}

static inline double logp_hd_u(enum State s, int idx, Intvl *intvl, int N)
{ Intvl I = intvl[idx];
  
  int nn_idx[2];
  find_nn_u(idx,s,intvl,N,nn_idx);
  int l_rel = nn_idx[0], r_rel = nn_idx[1];
  double logp_l, logp_r;
  
  { double logp_l_er = -INFINITY, logp_l_sf = -INFINITY, logp_l_sf_er = -INFINITY;
    int l = idx-1;
    if (l >= 0 && intvl[l].asgn == (char)s)
      logp_l_er = I.pe_o.b;
    if (l_rel != -1)
      { Intvl L = intvl[l_rel];
        logp_l_sf = logp_trans(L.e-1,I.b,L.cce,I.cb,L.cce);
      }
    cnt_t est_cnt = est_cov(I.b,idx,intvl,N,s,false);
    if (est_cnt >= I.cb)
      { double max_erate = 0.1;   // TODO: context
        logp_l_sf_er = log(p_errorin(OTHERS,max_erate,est_cnt,I.cb));
      }
    logp_l = MAX(MAX(logp_l_er,logp_l_sf),logp_l_sf_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
    fprintf(stderr,"  [%c L] logp(SF)=%lf, logp(ER)=%lf, logp(SF-ER)=%lf (-> %lf)\n",
                  stoc[s],logp_l_sf,logp_l_er,logp_l_sf_er,logp_l);
#endif
  }

  { double logp_r_er = -INFINITY, logp_r_sf = -INFINITY, logp_r_sf_er = -INFINITY;
    int r = idx+1;
    if (r < N && intvl[r].asgn == (char)s)
      logp_r_er = I.pe_o.e;
    if (r_rel != -1)
      { Intvl R = intvl[r_rel];
        logp_r_sf = logp_trans(I.e-1,R.b,I.ce,R.ccb,R.ccb);
      }
    cnt_t est_cnt = est_cov(I.e-1,idx,intvl,N,s,false);
    if (est_cnt >= I.ce)
      { double max_erate = 0.1;   // TODO: context
        logp_r_sf_er = log(p_errorin(OTHERS,max_erate,est_cnt,I.ce));
      }
    logp_r = MAX(MAX(logp_r_er,logp_r_sf),logp_r_sf_er);

#if defined(DEBUG_UNREL) && defined(DEBUG_PROB)
    fprintf(stderr,"  [%c R] logp(SF)=%lf, logp(ER)=%lf, logp(SF-ER)=%lf (-> %lf)\n",
                  stoc[s],logp_r_sf,logp_r_er,logp_r_sf_er,logp_r);
#endif
  }

  if (logp_l == -INFINITY && logp_r == -INFINITY)
    { logp_l = logp_poisson(I.cb,GLOBAL_COV[s]);
      logp_r = logp_poisson(I.ce,GLOBAL_COV[s]);
    }
  else if (logp_l == -INFINITY)
    logp_l = logp_r;
  else if (logp_r == -INFINITY)
    logp_r = logp_l;
  double logp = logp_l+logp_r;
  return logp;
}

static inline double logp_h_u(int idx, Intvl *intvl, int N)
{ return logp_hd_u(HAPLO,idx,intvl,N);
}

static inline double logp_d_u(int idx, Intvl *intvl, int N)
{ return logp_hd_u(DIPLO,idx,intvl,N);
}

static inline double calc_logp_u(enum State s, int idx, Intvl *intvl, int N)
{ if (s == ERROR) return logp_e_u(idx,intvl);
  else if (s == HAPLO) return logp_h_u(idx,intvl,N);
  else if (s == DIPLO) return logp_d_u(idx,intvl,N);
  else return logp_r_u(idx,intvl,N);
}

static void update_state(int idx, Intvl *intvl, int N
#ifdef DEBUG_UNREL
                         , int d
#endif
                         )
{ Intvl I = intvl[idx];

  if (MAX(I.cb,I.ce) >= GLOBAL_COV[REPEAT])
    { intvl[idx].asgn = REPEAT;
      return;
    }

#ifdef DEBUG_UNREL
      fprintf(stderr,"\n<%d> Intvl[%d] @ (%d,%d) %d -> %d: class=%c\n",
                     d,idx,I.b,I.e,I.cb,I.ce,(I.asgn == N_STATE) ? '-' : stoc[(int)I.asgn]);
#endif

  double logpmax = -INFINITY;
  int    smax = -1;
  for (int s = ERROR; s <= DIPLO; s++)
    { double logp = calc_logp_u(s,idx,intvl,N);
      if (logpmax < logp)
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

  if (I.asgn != smax)
    {
#ifdef DEBUG_UNREL
      fprintf(stderr,"state updated @ %d: %c -> %c\n",
                     idx,(I.asgn == N_STATE) ? '-' : stoc[(int)I.asgn],stoc[smax]);
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

void classify_unrel(Intvl *intvl, int N)
{ bool is_fixed[N];
  for (int i = 0; i < N; i++)
    is_fixed[i] = (intvl[i].is_rel && (intvl[i].asgn == HAPLO || intvl[i].asgn == DIPLO));

  Intvl_IC iord[N];
  for (int i = 0; i < N; i++)
    { iord[i].idx = i;
      iord[i].cnt = MIN(intvl[i].cb,intvl[i].ce);
    }
  qsort(iord,N,sizeof(Intvl_IC),compare_iic);

  for (int i = N-1; i >= 0; i--)
    if (!is_fixed[iord[i].idx])
      update_state(iord[i].idx,intvl,N
#ifdef DEBUG_UNREL
                   ,0
#endif
                   );

  for (int i = 0; i < N; i++)
    if (!is_fixed[iord[i].idx])
      update_state(iord[i].idx,intvl,N
#ifdef DEBUG_UNREL
                   ,1
#endif
                   );

#ifdef DEBUG_UNREL
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
  fprintf(stderr,"  Final: ");
  for (int i = 0; i < N; i++)
    for (int j = intvl[i].b; j < intvl[i].e; j++)
      fprintf(stderr,"%c",stoc[(int)intvl[i].asgn]);
  fprintf(stderr,"\n");
  fflush(stderr);
#endif
  return;
}

// TODO: debug; make it work
// void remove_slip(uint16 *profile, int plen, Seq_Ctx *ctx[N_WTYPE], char *crack)
// { int b, cnt, prev_d;
//   bool in_hp;   // TODO: extend to DS/TS?
//   char new_s;

//   b = 0;
//   cnt = 0;
//   prev_d = -1;
//   for (int i = 1; i < plen; i++)
//     { int d = (profile[i-1] >= profile[i]) ? DROP : GAIN;
//       if (d == DROP)
//         in_hp = (ctx[d][i-1][HP] < ctx[d][i][HP]) ? true : false;
//       else
//         in_hp = (ctx[d][i-1][HP] > ctx[d][i][HP]) ? true : false;
//       if (!in_hp || prev_d != d)
//         { cnt = 0;
//           continue;
//         }
//       if (crack[i-1] != crack[i])
//         { cnt++;
//           if (cnt == 1)
//             b = i;
//           else if (cnt >= 2)
//             { if (d == DROP)
//                 new_s = (MIN(crack[i-1],crack[i]) == ERROR) ? ERROR : crack[i-1];
//               else
//                 new_s = (MIN(crack[b-1],crack[b]) == ERROR) ? ERROR : crack[b];
// #ifdef DEBUG_SLIP
//               fprintf(stderr,"SLIP: [%d,%d) -> %c\n",b,i,stoc[(unsigned char)new_s]);
// #endif
//               for (int j = b; j < i; j++)
//                 crack[j] = new_s;
//             }
//         }
//       prev_d = d;
//     }

//   return;
// }
