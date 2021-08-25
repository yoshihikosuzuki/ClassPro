#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "ClassPro.h"

const int MAX_N_HC       = 5;   // Max # of bases in a single high-complexity error event
int       MIN_CNT_CHANGE = 2;   // Every count change at a wall must be > this
int       MAX_CNT_CHANGE = 5;   // Every count change > this becomes a wall candidate(H-cov?)

static const double THRES_LOGP_DIFF = log(1e-20);

const double pethres_init[N_ETYPE] = {0.001, 0.05};
const double pethres[N_ETYPE]      = {1e-10, 1e-5};

static inline double logp_diff(int cout_drop, int cin_drop, int cin_gain, int cout_gain)
{ int ndrop = cout_drop - cin_drop;
  int ngain = cout_gain - cin_gain;
  double _lambda = (double)MAX(cout_drop,cout_gain)/READ_LEN;
  return logp_skellam(ndrop-ngain,_lambda);
}

static void check_drop(int i, const uint16 *profile, int plen, Seq_Ctx *ctx[2], Error_Model *emodel, P_Error *perror, Error_Intvl *eintvl[N_ETYPE], const int eidx, const int K)
{ double max_p[N_ETYPE] = {-1., -1.};
  int    max_j[N_ETYPE] = {-1, -1};

  const double _pe = emodel[0].pe[1];
  double pe, ps, po, logp_d;
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
      // check p_diff (i.e. Skellam)?
      logp_d = logp_diff(profile[i-1],profile[i],profile[j-1],profile[j]);

      if (max_p[0] < ps)
        { max_p[0] = ps;
          max_j[0] = j;
        }
      if (logp_d > THRES_LOGP_DIFF && max_p[1] < po)
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
      fprintf(stderr,"    Po (= %lf * %lf * %lf) = %lf (LPd = %lf)",perror[i][OTHERS][DROP],perror[j][OTHERS][GAIN],pe,po,logp_d);
      if (logp_d > THRES_LOGP_DIFF && po > pethres[OTHERS])
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

          if (profile[j-1]>=profile[j])   // NOTE: not allowing tie counts
            { perror[j][SELF][GAIN] = 0.;
              perror[j][OTHERS][GAIN] = 1.;
            }
          else
            { perror[j][SELF][GAIN] = binom_test_g(profile[j-1],profile[j],emodel[t].pe[MIN(ctx[1][j][t],emodel[t].lmax)],0);
              perror[j][OTHERS][GAIN] = binom_test_g(profile[j]-profile[j-1],profile[j],emodel[t].pe[MIN(ctx[1][j][t],emodel[t].lmax)],0);
            }

          ps = perror[i][SELF][DROP] * perror[j][SELF][GAIN];
          po = perror[i][OTHERS][DROP] * perror[j][OTHERS][GAIN];
          logp_d = logp_diff(profile[i-1],profile[i],profile[j-1],profile[j]);

#ifdef DEBUG_ERROR
      fprintf(stderr,"  @ j = %d (t = %d, m = %d, n = %d): %d -> %d\n",j,t,m,n,profile[j-1],profile[j]);
      fprintf(stderr,"    Ps (= %lf * %lf) = %lf",perror[i][SELF][DROP],perror[j][SELF][GAIN],ps);
      if (ps > pethres[SELF])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
      fprintf(stderr,"    Po (= %lf * %lf) = %lf (LPd = %lf)",perror[i][OTHERS][DROP],perror[j][OTHERS][GAIN],po,logp_d);
      if (logp_d > THRES_LOGP_DIFF && po > pethres[OTHERS])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
#endif
        }

      if (max_p[0] < ps)
        { max_p[0] = ps;
          max_j[0] = j;
        }
      if (logp_d > THRES_LOGP_DIFF && max_p[1] < po)
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

static void check_gain(int i, const uint16 *profile, Seq_Ctx *ctx[2], Error_Model *emodel, P_Error *perror, Error_Intvl *eintvl[N_ETYPE], const int eidx, const int K)
{ double max_p[N_ETYPE] = {-1., -1.};
  int    max_j[N_ETYPE] = {-1, -1};

  const double _pe = emodel[0].pe[1];
  double pe, ps, po, logp_d;
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
      // check p_diff (i.e. Skellam)?
      logp_d = logp_diff(profile[j-1],profile[j],profile[i-1],profile[i]);
      // TODO: Multiply Pr(read start)^(# of count differences)?

      if (max_p[0] < ps)
        { max_p[0] = ps;
          max_j[0] = j;
        }
      if (logp_d > THRES_LOGP_DIFF && max_p[1] < po)
        { max_p[1] = po;
          max_j[1] = j;
        }

#ifdef DEBUG_ERROR
      fprintf(stderr,"  @ j = %d (m = %d, n = %d): %d -> %d\n",j,m,n,profile[j-1],profile[j]);
      fprintf(stderr,"    Ps (= %lf * %lf * %lf) = %lf",perror[j][SELF][DROP],perror[i][SELF][GAIN],pe,ps);
      if (ps > pethres[SELF])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
      fprintf(stderr,"    Po (= %lf * %lf * %lf) = %lf (LPd = %lf)",perror[j][OTHERS][DROP],perror[i][OTHERS][GAIN],pe,po,logp_d);
      if (logp_d > THRES_LOGP_DIFF && po > pethres[OTHERS])
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

          if (profile[j-1]<=profile[j])   // not allowing tie counts
            { perror[j][SELF][DROP] = 0.;
              perror[j][OTHERS][DROP] = 1.;
            }
          else
            { perror[j][SELF][DROP] = binom_test_g(profile[j],profile[j-1],emodel[t].pe[MIN(ctx[0][j][t],emodel[t].lmax)],0);
              perror[j][OTHERS][DROP] = binom_test_g(profile[j-1]-profile[j],profile[j-1],emodel[t].pe[MIN(ctx[0][j][t],emodel[t].lmax)],0);
            }

          ps = perror[j][SELF][DROP] * perror[i][SELF][GAIN];
          po = perror[j][OTHERS][DROP] * perror[i][OTHERS][GAIN];
          logp_d = logp_diff(profile[j-1],profile[j],profile[i-1],profile[i]);

#ifdef DEBUG_ERROR
      fprintf(stderr,"  @ j = %d (t = %d, m = %d, n = %d): %d -> %d\n",j,t,m,n,profile[j-1],profile[j]);
      fprintf(stderr,"    Ps (= %lf * %lf) = %lf",perror[j][SELF][DROP],perror[i][SELF][GAIN],ps);
      if (ps > pethres[SELF])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
      fprintf(stderr,"    Po (= %lf * %lf) = %lf (LPd = %lf)",perror[j][OTHERS][DROP],perror[i][OTHERS][GAIN],po,logp_d);
      if (logp_d > THRES_LOGP_DIFF && po > pethres[OTHERS])
        fprintf(stderr," ***");
      fprintf(stderr,"\n");
#endif
        }

      if (max_p[0] < ps)
        { max_p[0] = ps;
          max_j[0] = j;
        }
      if (logp_d > THRES_LOGP_DIFF && max_p[1] < po)
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

static void correct_wall_cnt(int beg, int end, int *ci, int *cj, const uint16 *profile, Seq_Ctx *ctx[N_WTYPE], const int K)
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

void find_wall(const uint16 *profile, int plen, Seq_Ctx *ctx[N_WTYPE], Error_Model *emodel, P_Error *perror, P_Error *cerror, Error_Intvl *eintvl[N_ETYPE], Intvl *intvl, Rel_Intvl *rintvl, int *wall, char *asgn, const int K, int *_N, int *_M)
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

      asgn[i] = 1;   // TODO: check if this is OK? is_wall <==> asgn=1

      if (perror[i][SELF][w] == 0.)
        { //perror[i][SELF][w] = (cout <= CMAX) ? emodel[maxt].pe_bt[maxl][CIDX(cout,cin)] : binom_test_g(cin,cout,maxpe,0);
          perror[i][SELF][w] = binom_test_g(cin,cout,maxpe,0);
        }
      if (perror[i][OTHERS][w] == 0.)
        { //perror[i][OTHERS][w] = (cout <= CMAX) ? emodel[maxt].pe_bt[maxl][CIDX(cout,cout-cin)] : binom_test_g(cout-cin,cout,maxpe,0);
          perror[i][OTHERS][w] = binom_test_g(cout-cin,cout,maxpe,0);
        }

#ifdef DEBUG_ERROR
      const char _type = (profile[i-1] == profile[i]) ? '=' : ((profile[i-1] > profile[i]) ? '>' : '<');
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
        { // asgn[i] = 1;

          if (w == DROP)
            { for (int n = 0; n <= MAX_N_HC; n++)
                { int j = i-1+K+n;
                  // if (j >= plen || profile[j-1]>profile[j])
                  if (j >= plen || profile[j-1] >= profile[j])   // NOTE: experimental code to skip tie
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
                  // if (j <= 0 || profile[j-1]<profile[j])
                  if (j <= 0 || profile[j-1] <= profile[j])   // NOTE: experimental code to skip tie
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
  // e = SELF;
  // for (int i = 0; i < eidx; i++)
  //   if (eintvl[e][i].pe >= pethres[e])
  //     asgn[eintvl[e][i].i] = asgn[eintvl[e][i].j] = 1;
  
  // Concatenate adjacent error intervals
  bool is_error[plen];
  for (int i = 0; i < plen; i++)
    is_error[i] = false;
  for (int i = 0; i < eidx; i++)
    if (eintvl[SELF][i].pe >= pethres[SELF])
      { //fprintf(stderr,"Eintvl[%d] = (%d,%d)\n",i,eintvl[SELF][i].i,eintvl[SELF][i].j);
        for (int j = eintvl[SELF][i].i; j < eintvl[SELF][i].j; j++)
          is_error[j] = true;
      }
  for (int i = 1; i < plen; i++)
    if (is_error[i-1] != is_error[i])
      asgn[i] = 1;
    else if (is_error[i-1] && is_error[i])
      asgn[i] = 0;

  // Store logp error in self considering the pair, for each interval
  double pe_i[plen+1], pe_j[plen+1];
  for (int i = 0; i < plen+1; i++)
    { pe_i[i] = pe_j[i] = -INFINITY;
    }
  for (int i = 0; i < eidx; i++)
    { // fprintf(stderr,"*** %d (%d,%d)=%lf\n",i,eintvl[SELF][i].i,eintvl[SELF][i].j,eintvl[SELF][i].pe);
      pe_i[eintvl[SELF][i].i] = MAX(pe_i[eintvl[SELF][i].i], log(eintvl[SELF][i].pe));
      pe_j[eintvl[SELF][i].j] = MAX(pe_i[eintvl[SELF][i].j], log(eintvl[SELF][i].pe));
    }

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
          fprintf(stderr,"(%d,%d) [%d, %d] Psk=%lf",beg,end,ci,cj,exp(logp_sk));
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
      intvl[iidx].logpe_i = pe_i[beg];
      intvl[iidx].logpe_j = pe_j[end];
      // fprintf(stderr,"PeS=(%lf,%lf)\n",pe_i[beg],pe_j[end]);
      intvl[iidx].is_rel = is_rel;
      intvl[iidx].is_err = is_err;
      intvl[iidx].asgn = (is_err) ? ERROR : N_STATE;
      iidx++;
    }

  *_N = N;
  *_M = ridx;
  return;
}
