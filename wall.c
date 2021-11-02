#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "ClassPro.h"
#include <gsl/gsl_multifit.h>

static void polynomialfit(int N, int degree, double *data_x, double *data_y, double *coef)
{ gsl_matrix *X, *cov;
  gsl_vector *y, *c;
  double chisq;

  X = gsl_matrix_alloc(N, degree);
  y = gsl_vector_alloc(N);
  c = gsl_vector_alloc(degree);
  cov = gsl_matrix_alloc(degree, degree);

  for (int i = 0; i < N; i++)
    { for (int j = 0; j < degree; j++)
        gsl_matrix_set(X,i,j,pow(data_x[i],j));
      gsl_vector_set(y,i,data_y[i]);
    }

  { gsl_multifit_linear_workspace *work;
    work = gsl_multifit_linear_alloc(N,degree);
    gsl_multifit_linear(X,y,c,cov,&chisq,work);
    gsl_multifit_linear_free(work);
  }

  for (int i = 0; i < degree; i++)
    coef[i] = gsl_vector_get(c,i);

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return;
}

typedef struct
  { float  all;
    float  ins;
    float  op[9];
  } E_Rates;

typedef struct
  { float all;
    float op[6];
  } M_Rates;

static void load_himodel(Error_Model *emodel, const char *name)
{ FILE *efile = fopen(name,"r");
  int   i, kmer;
  fread(&kmer,sizeof(int),1,efile);
  const int krange = kmer/2 - 6;

  E_Rates HepTab[0x4000];
  fread(HepTab,sizeof(E_Rates),0x4000,efile);

  M_Rates *mics[N_CTYPE];
  int KMAX[N_CTYPE];
  double coef[3];
  double x[5] = {1, 2, 3, 4, 5}, y[5];
  int n[5];
  y[0] = 0.002;
  for (int t = HP; t <= TS; t++)
    { int ulen = t+1;
      int N = (1 << (2*ulen));
      mics[t] = ((M_Rates *)Malloc(sizeof(M_Rates)*N*krange,"Table"));
      fread(mics[t],sizeof(M_Rates),krange*N,efile);
      mics[t] -= 2*ulen;
      KMAX[t] = 2*ulen+krange;

      if (VERBOSE)
        { for (i = 0; i < N; i++)
            { printf("ulen=%d, i=%d:\n",ulen,i);
              for (int j = 2*ulen; j < KMAX[t]; j++)
                printf("  [%d] %5.2f",j,100.*mics[t][krange*i+j].all);
              printf("\n");
            }
        }

      for (int j = 2; j <= 5; j++)
        { y[j-1] = 0.;
          n[j-1] = 0;
        }
      for (int j = 2; j <= 5; j++)
        { for (i = 0; i < N; i++)
            { double p = mics[t][krange*i+j*ulen].all;
              if (p > 0.)
                { y[j-1] += p;
                  n[j-1]++;
                }
            }
          y[j-1] /= n[j-1];
        }
      polynomialfit(5,3,x,y,coef);
      for (int l = 1; l <= emodel[t].lmax; l++)
        emodel[t].pe[l] = coef[0]+coef[1]*l+coef[2]*l*l;

      if (VERBOSE)
        { for (int j = 2*ulen; j < KMAX[t]; j++)
            { int _j = j/ulen;
              printf("  [%d] %5.2f",j,100.*(coef[0]+coef[1]*_j+coef[2]*_j*_j));
            }
          printf(" (%f + %f x + %f x^2)\n",coef[0],coef[1],coef[2]);
        }
    }

  return;
}

static uint8  CMAX;
static double HC_ERATE;

static Error_Model *load_emodel(const char *name)
{ Error_Model *emodel = Malloc(sizeof(Error_Model)*N_CTYPE,"Allocating error model");
  for (int t = HP; t <= TS; t++)
    { emodel[t].lmax = (uint8)(MAX_N_LC/(t+1));
      emodel[t].pe = Malloc(sizeof(double)*(emodel[t].lmax+1),"Allocating pe");
      emodel[t].pe[0] = 0.;
      emodel[t].cthres = Malloc(sizeof(uint8***)*(emodel[t].lmax+1),"Allocating cthres");
      for (int l = 1; l <= emodel[t].lmax; l++)
        { emodel[t].cthres[l] = Malloc(sizeof(uint8**)*CMAX,"Allocating cthres l");
          for (int c = 1; c < CMAX; c++)
            { emodel[t].cthres[l][c] = Malloc(sizeof(uint8*)*N_THRES,"Allocating cthres pe");
              for (int s = INIT; s <= FINAL; s++)
                emodel[t].cthres[l][c][s] = Malloc(sizeof(uint8)*N_ETYPE,"Allocating cthres thresT");
            }
        }
    }

  if (name == NULL)
    { if (VERBOSE)
        fprintf(stderr,"Error model not specified. Use the default error model.\n");
      for (int t = HP; t <= TS; t++)
        for (int l = 1; l <= emodel[t].lmax; l++)
          emodel[t].pe[l] = 0.002 * l * l + 0.002;
    }
  else
    load_himodel(emodel,name);

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

Error_Model *calc_init_thres(const char *name)
{ Error_Model *emodel;
  uint8        cout, cin;
  uint8        ct[N_ETYPE];
  double       psum;
  bool         is_found[N_THRES][N_ETYPE];

  if (GLOBAL_COV[REPEAT] > 255)
    { fprintf(stderr,"Too high REPEAT coverage (%d) > 255\n",GLOBAL_COV[REPEAT]);
      exit(1);
    }
  CMAX = GLOBAL_COV[REPEAT];
  emodel = load_emodel(name);
  HC_ERATE = emodel[HP].pe[1];

#ifdef DEBUG_EMODEL
  fprintf(stderr,"Thresholds for initial wall filtering (cmax = %d):\n",CMAX);
  fprintf(stderr,"          cout       :");
  for (cout = 1; cout < CMAX; cout++)
    fprintf(stderr," %3d",cout);
  fprintf(stderr,"\n  ( t, l, pe)\n");
#endif

  for (int t = HP; t <= TS; t++)
    for (int l = 1; l <= emodel[t].lmax; l++)
      { double pe = emodel[t].pe[l];
        double lpe = log(pe);
        double l1mpe = log(1-pe);

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

Wall_Arg *alloc_wall_arg(int rlen_max)
{ Wall_Arg *arg = Malloc(sizeof(Wall_Arg),"Allocating wall arg");
  arg->wall     = Malloc((rlen_max+1)*sizeof(char),"Wall array");
  arg->eintvl   = Malloc(rlen_max*sizeof(Error_Intvl),"Error(S) intvl array");
  arg->ointvl   = Malloc(rlen_max*sizeof(Error_Intvl),"Error(O) intvl array");
  arg->perror   = Malloc((rlen_max+1)*sizeof(P_Error),"Error prob array");
  return arg;
}

void free_wall_arg(Wall_Arg *arg)
{ free(arg->wall);
  free(arg->eintvl);
  free(arg->ointvl);
  free(arg->perror);
  free(arg);
  return;
}

static const char MASK_WALL_BY[N_ETYPE]   = { 0x01, 0x10 };
static const char MASK_PAIRED_BY[N_ETYPE] = { 0x02, 0x20 };
static const char MASK_WALL_MULT          = 0x04;
static const char MASK_PAIRED_MULT        = 0x40;
static const char MASK_WALL               = 0x08;
static const char MASK_ERROR              = 0x80;

static inline bool is_wall_by(enum Etype e, char *wall, int i)
  { return (wall[i] & MASK_WALL_BY[e]) ? true : false; }

static inline void set_wall_by(enum Etype e, char *wall, int i)
  { wall[i] |= MASK_WALL_BY[e]; return; }

static inline void unset_wall_by(enum Etype e, char *wall, int i)
  { wall[i] &= ~(MASK_WALL_BY[e]); return; }

static inline bool is_paired(enum Etype e, char *wall, int i)
  { return (wall[i] & MASK_PAIRED_BY[e]) ? true : false; }

static inline void set_paired(enum Etype e, char *wall, int i)
  { wall[i] |= MASK_PAIRED_BY[e]; return; }

static inline bool is_wall_by_mult(char *wall, int i)
  { return (wall[i] & MASK_WALL_MULT) ? true : false; }

static inline void set_wall_by_mult(char *wall, int i)
  { wall[i] |= MASK_WALL_MULT; return; }

static inline bool is_paired_mult(char *wall, int i)
  { return (wall[i] & MASK_PAIRED_MULT) ? true : false; }

static inline void set_paired_mult(char *wall, int i)
  { wall[i] |= MASK_PAIRED_MULT; return; }

static inline bool is_wall(char *wall, int i)
  { return (wall[i] & MASK_WALL) ? true : false; }

static inline void set_wall(char *wall, int i)
  { wall[i] |= MASK_WALL; return; }

static inline bool is_error(char *wall, int i)
  { return (wall[i] & MASK_ERROR) ? true : false; }

static inline void set_error(char *wall, int i)
  { wall[i] |= MASK_ERROR; return; }

static inline void update_perror(P_Error *perror, int i, enum Etype e, enum Wtype w,
                                 cnt_t cout, cnt_t cin, double erate)
{ if (perror[i][e][w] == -INFINITY)
    perror[i][e][w] = p_errorin(e,erate,cout,cin);
  return;
}

static inline double logp_diff_pair(pos_t i, pos_t j, const cnt_t *profile)
{ int n_drop = (int)profile[i-1]-profile[i];
  int n_gain = (int)profile[j]-profile[j-1];
  cnt_t cov  = MAX(profile[i-1],profile[j]);
  return logp_trans(i,j,n_drop,n_gain,cov);
}

static inline bool cthres_ng(enum Etype e, uint8 cin, uint8 ct)
{ if (e == SELF)
    return (cin >= ct) ? true : false;
  else
    return (cin < ct) ? true : false;
}

static bool find_gain(int i, cnt_t cout, cnt_t cin, enum Etype e, enum Ctype t, int l, double erate,
                      P_Error *perror, Error_Intvl *intvl, int *idx,
                      cnt_t *profile, int plen, Seq_Ctx *ctx[N_WTYPE], Error_Model *emodel, int K)
{ int    m, n, j, max_j;
  double pe, max_pe;
  cnt_t  cout_j, cin_j;
  
  const int ipk  = i+K-1;
  const int ulen = t+1;

  max_j = -1;
  max_pe = -INFINITY;

  // Low-complexity error
  m = ulen*l;
  n = 0;
  while (true)
    { int idx = i+ulen*(n+1);
      if (idx >= plen || ctx[DROP][idx][t] != m+n+1)
        break;
      n++;
    }
  j = ipk+n-m;
  if (j <= i)
    return false;
  if (j >= plen)
    { j = plen;
      pe = perror[i][e][DROP] * perror[i][e][DROP];
    }
  else
    { cin_j  = profile[j-1];
      cout_j = profile[j];
      pe = -INFINITY;
      if (cin_j <= cout_j
          && !(cout_j < CMAX && cthres_ng(e,cin_j,emodel[t].cthres[l][cout_j][FINAL][e]))
          && (e == SELF || logp_diff_pair(i,j,profile) >= THRES_DIFF_EO))
        { update_perror(perror,j,e,GAIN,cout_j,cin_j,erate);
          pe = perror[i][e][DROP]*perror[j][e][GAIN];
        }
    }
  if (max_pe < pe)
    { max_j  = j;
      max_pe = pe;
    }

  // High-complexity errors
  m = 0;
  for (n = 0; n <= MAX_N_HC; n++)
    { j = ipk+n-m;
      if (j >= plen)
        break;
      cin_j  = profile[j-1];
      cout_j = profile[j];
      if (!(cin_j <= cout_j))
        continue;

      if ((cout < CMAX && cthres_ng(e,cin,emodel[HP].cthres[1][cout][FINAL][e]))
          || (cout_j < CMAX && cthres_ng(e,cin_j,emodel[HP].cthres[1][cout_j][FINAL][e])))
        continue;
      if (e == OTHERS && logp_diff_pair(i,j,profile) < THRES_DIFF_EO)
        continue;
      
      // TODO: improve this heavy computation with small gain...
      double pe_i = p_errorin(e,HC_ERATE,cout,cin);
      double pe_j = p_errorin(e,HC_ERATE,cout_j,cin_j);
      pe = pe_i * pe_j;
      if (max_pe < pe)
        { max_j  = j;
          max_pe = pe;
        }
    }

  if (max_j == -1)
    return false;

  intvl[*idx].b  = i;
  intvl[*idx].e  = max_j;
  intvl[*idx].pe = max_pe;

// #ifdef DEBUG_WALL
//   fprintf(stderr,"  @ max_j = %d (d=%d): %d -> %d (pe=%lf)\n",
//                  max_j,max_j-i,profile[max_j-1],profile[max_j],max_pe);
// #endif

  return true;
}

static bool find_drop(int i, cnt_t cout, cnt_t cin, enum Etype e, enum Ctype t, int l, double erate,
                      P_Error *perror, Error_Intvl *intvl, int *idx,
                      cnt_t *profile, Seq_Ctx *ctx[N_WTYPE], Error_Model *emodel, int K)
{ int    m, n, j, max_j;
  double pe, max_pe;
  cnt_t  cout_j, cin_j;
  
  const int imk  = i-K+1;
  const int ulen = t+1;

  max_j = -1;
  max_pe = -INFINITY;

  // Low-complexity error
  m = ulen*l;
  n = 0;
  while (true)
    { int idx = i-ulen*(n+1);
      if (idx < 0 || ctx[GAIN][idx][t] != m+n+1)
        break;
      n++;
    }
  j = imk-n+m;
  if (j >= i)
    return false;
  if (j <= 0)
    { j = 0;
      pe = perror[i][e][GAIN] * perror[i][e][GAIN];
    }
  else
    { cout_j = profile[j-1];
      cin_j  = profile[j];
      pe = -INFINITY;
      if (cin_j <= cout_j
          && !(cout_j < CMAX && cthres_ng(e,cin_j,emodel[t].cthres[l][cout_j][FINAL][e]))
          && (e == SELF || logp_diff_pair(j,i,profile) >= THRES_DIFF_EO))
        { update_perror(perror,j,e,DROP,cout_j,cin_j,erate);
          pe = perror[j][e][DROP]*perror[i][e][GAIN];
        }
    }
  if (max_pe < pe)
    { max_j  = j;
      max_pe = pe;
    }

  // High-complexity errors
  m = 0;
  for (n = 0; n <= MAX_N_HC; n++)
    { j = imk-n+m;
      if (j <= 0)
        break;
      cout_j = profile[j-1];
      cin_j  = profile[j];
      if (!(cin_j <= cout_j))
        continue;

      if ((cout < CMAX && cthres_ng(e,cin,emodel[HP].cthres[1][cout][FINAL][e]))
          || (cout_j < CMAX && cthres_ng(e,cin_j,emodel[HP].cthres[1][cout_j][FINAL][e])))
        continue;
      if (e == OTHERS && logp_diff_pair(j,i,profile) < THRES_DIFF_EO)
        continue;
      
      // TODO: improve this heavy computation with small gain...
      double pe_i = p_errorin(e,HC_ERATE,cout,cin);
      double pe_j = p_errorin(e,HC_ERATE,cout_j,cin_j);
      pe = pe_i * pe_j;
      if (max_pe < pe)
        { max_j  = j;
          max_pe = pe;
        }
    }

  if (max_j == -1)
    return false;

  intvl[*idx].b  = max_j;
  intvl[*idx].e  = i;
  intvl[*idx].pe = max_pe;

// #ifdef DEBUG_WALL
//   fprintf(stderr,"  @ max_j = %d (d=%d): %d -> %d (pe=%lf)\n",
//                  max_j,i-max_j,profile[max_j-1],profile[max_j],max_pe);
// #endif

  return true;
}

static inline bool find_pair(int i, cnt_t cout, cnt_t cin,
                             enum Etype e, enum Wtype w, enum Ctype t, int l, double erate,
                             P_Error *perror, Error_Intvl *intvl, int *idx,
                             cnt_t *profile, int plen, Seq_Ctx *ctx[N_WTYPE], Error_Model *emodel, int K)
{ if (w == DROP)
    return find_gain(i,cout,cin,e,t,l,erate,perror,intvl,idx,profile,plen,ctx,emodel,K);
  else
    return find_drop(i,cout,cin,e,t,l,erate,perror,intvl,idx,profile,ctx,emodel,K);
}

static int compare_eintvl(const void *a, const void *b)
{ if (((Error_Intvl*)a)->b == ((Error_Intvl*)b)->b)
    { if (((Error_Intvl*)a)->e == ((Error_Intvl*)b)->e)
        return ((Error_Intvl*)b)->pe - ((Error_Intvl*)a)->pe;
      else
        return ((Error_Intvl*)a)->e - ((Error_Intvl*)b)->e;
    }
  else
    return ((Error_Intvl*)a)->b - ((Error_Intvl*)b)->b;
}

static int bs_eintvl(Error_Intvl *intvl, int l, int r, pos_t b, pos_t e) {
  if (l > r)
     return -1;
  int m = (l+r)/2;
  if (intvl[m].b == b)
    { if (intvl[m].e == e)
        return m;
      else if (e > intvl[m].e)
        return bs_eintvl(intvl,m+1,r,b,e);
      else
        return bs_eintvl(intvl,l,m-1,b,e);
    }
  else if (b > intvl[m].b)
    return bs_eintvl(intvl,m+1,r,b,e);
  else
    return bs_eintvl(intvl,l,m-1,b,e);
}

static inline int remove_duplicates(Error_Intvl *intvl, int N)
{ qsort(intvl,N,sizeof(Error_Intvl),compare_eintvl);
  if (N >= 2)
    { int i = 1;
      while (i < N)
        { if (intvl[i-1].b == intvl[i].b && intvl[i-1].e == intvl[i].e)
            break;
          i++;
        }
      for (int j = i+1; j < N; j++)
        { if (!(intvl[i-1].b == intvl[j].b && intvl[i-1].e == intvl[j].e))
            { intvl[i].b  = intvl[j].b;
              intvl[i].e  = intvl[j].e;
              intvl[i].pe = intvl[j].pe;
              i++;
            }
        }
      N = i;
    }
  return N;
}

int find_wall(Wall_Arg *arg, Intvl *intvl, cnt_t *profile, int plen,
              Seq_Ctx *ctx[N_WTYPE], Error_Model *emodel, int K)
{ char        *wall      = arg->wall;
  Error_Intvl *eintvl    = arg->eintvl;   // TODO: merge eintvl and ointvl?
  Error_Intvl *ointvl    = arg->ointvl;
  P_Error     *perror    = arg->perror;

#ifdef DEBUG_WALL
  fprintf(stderr,"\n");
#endif

  for (int i = 0; i < plen; i++)
    { wall[i] = 0;
      for (int e = SELF; e <= OTHERS; e++)
        for (int w = DROP; w <= GAIN; w++)
          perror[i][e][w] = -INFINITY;
    }

  uint8 ct[N_THRES];
  int eidx = 0, oidx = 0;
  for (int i = 1; i < plen; i++)
    { // Determine wall type and set cout/cin
      cnt_t cim1 = profile[i-1];
      cnt_t ci   = profile[i];
      if (MIN(cim1,ci) >= GLOBAL_COV[REPEAT])
        continue;

      uint16 cng = abs((int)cim1-ci);
      if (cng < MIN_CNT_CHANGE)
        continue;

      enum Wtype wtype;
      cnt_t      cin, cout;
      if (cim1 > ci)
        { wtype = DROP;
          cin   = ci;
          cout  = cim1;
        }
      else
        { wtype = GAIN;
          cin   = cim1;
          cout  = ci;
        }

      // Determine low-complexity sequence context type
      int    maxt = -1, maxl = -1;
      double pe, maxpe = -INFINITY;
      for (int t = HP; t <= TS; t++)
        { int l = MIN(ctx[wtype][i][t],emodel[t].lmax);
          pe = emodel[t].pe[l];
          if (maxpe < pe)
            { maxpe = pe;
              maxt  = t;
              maxl  = l;
            }
        }

      for (int e = SELF; e <= OTHERS; e++)
        { if (is_paired(e, wall, i))
            continue;
          
          // Load precomputed count thresholds if available
          if (cout < CMAX)
            { for (int s = INIT; s <= FINAL; s++)
                ct[s] = emodel[maxt].cthres[maxl][cout][s][e];   // TODO: change to cthres[...][e][s] ?
              if (!(cng > MAX_CNT_CHANGE || cin < MAX(ct[INIT],3)))
                continue;
            }
          
          // Check if the position is wall
          if (e == SELF)
            { if (cout < CMAX && cin >= ct[FINAL])   // Cannot be E-intvl
                continue;
              update_perror(perror,i,e,wtype,cout,cin,maxpe);
              if (perror[i][e][wtype] < PE_THRES[FINAL][e])
                continue;
              if (find_pair(i,cout,cin,e,wtype,maxt,maxl,maxpe,perror,eintvl,&eidx,profile,plen,ctx,emodel,K))
                { Error_Intvl I = eintvl[eidx];
                  if (I.pe >= PE_THRES[FINAL][e])
                    { set_wall_by(e,wall,I.b);
                      set_wall_by(e,wall,I.e);
                      set_paired(e,wall,I.b);
                      set_paired(e,wall,I.e);
                      eidx++;
                    }
                }
              // TODO: need this?
              // else
              //   set_wall_by(e,wall,i);
            }
          else   // OTHERS
            { if (cng >= GLOBAL_COV[HAPLO] || (cout < CMAX && cin < ct[FINAL]))
                { set_wall_by(OTHERS,wall,i);
                  continue;
                }
              update_perror(perror,i,e,wtype,cout,cin,maxpe);
              if (perror[i][e][wtype] < PE_THRES[FINAL][e])   // Never paired
                { set_wall_by(OTHERS,wall,i);
                  continue;
                }
              if (find_pair(i,cout,cin,e,wtype,maxt,maxl,maxpe,perror,ointvl,&oidx,profile,plen,ctx,emodel,K))
                { Error_Intvl I = ointvl[oidx];
                  if (I.pe >= PE_THRES[FINAL][e])
                    { set_paired(e,wall,I.b);
                      set_paired(e,wall,I.e);
                      oidx++;
                      continue;
                    }
                }
              set_wall_by(e,wall,i);
            }
        }

// #ifdef DEBUG_WALL
//       if (is_wall_by(SELF,wall,i) || is_wall_by(OTHERS,wall,i))
//         { const char _type = (cim1 == ci) ? '=' : ((wtype == DROP) ? '>' : '<');
//           fprintf(stderr,"@ %d -> %d: %d -> %d (%c) ctx=",i-1,i,cim1,ci,_type);
//           for (int w = DROP; w <= GAIN; w++)
//             { fprintf(stderr,"%s(",(w == DROP) ? "L" : " R");
//               for (int t = HP; t <= TS; t++)
//                 fprintf(stderr,"%d%c",ctx[w][i][t],(t < TS) ? ',' : ')');
//             }
//           fprintf(stderr," erate=%.3lf, (ps, po) = (%lf, %lf)\n",
//                          maxpe,perror[i][SELF][wtype],perror[i][OTHERS][wtype]);
//         }
// #endif
    }   // Loop for each position

  int NS = eidx, NO = oidx;

#ifdef DEBUG_WALL
  int NS_init = NS, NO_init = NO;
  // fprintf(stderr,"E-intvls (init=%d):\n",NS);
  // for (int i = 0; i < NS; i++)
  //   fprintf(stderr,"    (%d, %d) pe=%lf\n",eintvl[i].b,eintvl[i].e,eintvl[i].pe);
  // fprintf(stderr,"O-intvls (init=%d):\n",NO);
  // for (int i = 0; i < NO; i++)
  //   fprintf(stderr,"    (%d, %d) pe=%lf\n",ointvl[i].b,ointvl[i].e,ointvl[i].pe);
#endif

  // From walls by non-O, exclude positions explained by O or within E-intvls
  for (int i = 0; i < NO; i++)
    { Error_Intvl I = ointvl[i];
      unset_wall_by(OTHERS,wall,I.b);
      unset_wall_by(OTHERS,wall,I.e);
    }
  for (int i = 0; i < NS; i++)
    { Error_Intvl I = eintvl[i];
      for (pos_t j = I.b+1; j < I.e; j++)
        unset_wall_by(OTHERS,wall,j);
    }

  // Sort by position and remove duplicates
  NS = remove_duplicates(eintvl,eidx);
  NO = remove_duplicates(ointvl,oidx);

#ifdef DEBUG_WALL
  fprintf(stderr,"E-intvls (init=%d, uniq=%d):\n",NS_init,NS);
  for (int i = 0; i < NS; i++)
    fprintf(stderr,"    (%d, %d) pe=%lf\n",eintvl[i].b,eintvl[i].e,eintvl[i].pe);
  fprintf(stderr,"O-intvls (init=%d, uniq=%d):\n",NO_init,NO);
  for (int i = 0; i < NO; i++)
    fprintf(stderr,"    (%d, %d) pe=%lf\n",ointvl[i].b,ointvl[i].e,ointvl[i].pe);
#endif

#ifdef DEBUG_WALL
  fprintf(stderr,"S-walls (init):");
  for (int i = 1; i < plen; i++)
    if (is_wall_by(SELF,wall,i))
      fprintf(stderr," %d",i);
  fprintf(stderr,"\n");
  fprintf(stderr,"O-walls (init):");
  for (int i = 1; i < plen; i++)
    if (is_wall_by(OTHERS,wall,i))
      fprintf(stderr," %d",i);
  fprintf(stderr,"\n");
#endif

  // Find E-intvls by multiple errors and boundary E-intvls
  int midx = NS;
  int b, e;
  double pe, pe_i, pe_j;
  for (int i = 1; i < plen; i++)
    { if (!(is_wall_by(OTHERS,wall,i) && !is_wall_by(SELF,wall,i)))
        continue;
      if (is_paired_mult(wall,i))
        continue;
      for (int w = DROP; w <= GAIN; w++)
        { if ((pe_i = perror[i][SELF][w]) < PE_THRES[FINAL][SELF])
            continue;
          if (w == DROP)
            { for (int j = i+1; j < MIN(i+200,plen+1); j++)
                { if (j == plen)   // boundary E-intvl
                    { if ((pe = pe_i * pe_i) < PE_THRES[FINAL][SELF])
                        continue;
                      b = i;
                      e = plen;
                      eintvl[midx].b = b;
                      eintvl[midx].e = e;
                      eintvl[midx].pe = pe;
                      set_paired_mult(wall,i);
                      midx++;
#ifdef DEBUG
                      if (midx >= plen)
                        { fprintf(stderr,"# E-intvls >= plen\n");
                          exit(1);
                        }
#endif
                    }
                  if (!is_wall_by(SELF,wall,j) && !is_wall_by(OTHERS,wall,j))
                    continue;
                  b = i;
                  e = j;
                  if (bs_eintvl(eintvl,0,NS,b,e) == -1)
                    { pe_j = perror[j][SELF][GAIN];
                      if ((pe = pe_i * pe_j) >= PE_THRES[FINAL][SELF])
                        { eintvl[midx].b = b;
                          eintvl[midx].e = e;
                          eintvl[midx].pe = pe;
                          set_paired_mult(wall,i);
                          set_paired_mult(wall,j);
                          midx++;
#ifdef DEBUG
                          if (midx >= plen)
                            { fprintf(stderr,"# E-intvls >= plen\n");
                              exit(1);
                            }
#endif
                        }
                    }
                  if (is_wall_by(OTHERS,wall,j))
                    break;
                }
            }
          else   // GAIN
            { for (int j = i-1; j >= MAX(i-200,0); j--)
                { if (j == 0)   // boundary E-intvl
                    { if ((pe = pe_i * pe_i) < PE_THRES[FINAL][SELF])
                        continue;
                      b = 0;
                      e = i;
                      eintvl[midx].b = b;
                      eintvl[midx].e = e;
                      eintvl[midx].pe = pe;
                      set_paired_mult(wall,i);
                      midx++;
#ifdef DEBUG
                      if (midx >= plen)
                        { fprintf(stderr,"# E-intvls >= plen\n");
                          exit(1);
                        }
#endif
                    }
                  if (!is_wall_by(SELF,wall,j) && !is_wall_by(OTHERS,wall,j))
                    continue;
                  b = j;
                  e = i;
                  if (bs_eintvl(eintvl,0,NS,b,e) == -1)
                    { pe_j = perror[j][SELF][DROP];
                      if ((pe = pe_i * pe_j) >= PE_THRES[FINAL][SELF])
                        { eintvl[midx].b = b;
                          eintvl[midx].e = e;
                          eintvl[midx].pe = pe;
                          set_paired_mult(wall,i);
                          set_paired_mult(wall,j);
                          midx++;
#ifdef DEBUG
                          if (midx >= plen)
                            { fprintf(stderr,"# E-intvls >= plen\n");
                              exit(1);
                            }
#endif
                        }
                    }
                  if (is_wall_by(OTHERS,wall,j))
                    break;
                }
            }
        }
    }
#ifdef DEBUG_WALL
  fprintf(stderr,"E-intvls (mult=%d):\n",midx-NS);
  for (int i = NS; i < midx; i++)
    fprintf(stderr,"    (%d, %d) pe=%lf\n",eintvl[i].b,eintvl[i].e,eintvl[i].pe);
#endif

  for (int i = NS; i < midx; i++)
    { Error_Intvl I = eintvl[i];
      for (pos_t j = I.b+1; j < I.e; j++)
        unset_wall_by(OTHERS,wall,j);
    }
  if (NS < midx)
    { NS = midx;
      qsort(eintvl,NS,sizeof(Error_Intvl),compare_eintvl);
    }

  // Merge overlapping/contained E-intvls (original E-intvls are kept)
  { int i = 0, j;
    pos_t max_e;
    double max_pe;
    while (i < NS-1)
      { max_e = eintvl[i].e;
        max_pe = eintvl[i].pe;
        j = i;
        while(j < NS-1)
          { if (eintvl[j+1].b <= eintvl[j].e)
              { max_e = MAX(max_e,eintvl[j+1].e);
                max_pe = MAX(max_pe,eintvl[j+1].pe);
                j++;
              }
            else
              break;
          }
        if (i < j)
          { eintvl[NS].b = eintvl[i].b;
            eintvl[NS].e = max_e;
            eintvl[NS].pe = max_pe;
            NS++;
#ifdef DEBUG
            if (NS >= plen)
              { fprintf(stderr,"# E-intvls >= plen\n");
                exit(1);
              }
#endif
          }
        i = j+1;
      }
  }
  qsort(eintvl,NS,sizeof(Error_Intvl),compare_eintvl);

// #ifdef DEBUG_WALL
//   fprintf(stderr,"E-intvls (all=%d):\n",NS);
//   for (int i = 0; i < NS; i++)
//     fprintf(stderr,"    (%d, %d) pe=%lf\n",eintvl[i].b,eintvl[i].e,eintvl[i].pe);
// #endif

  for (int i = 0; i < NS; i++)
    for (pos_t j = eintvl[i].b; j < eintvl[i].e; j++)
      set_error(wall,j);

  // Determine intervals
  int _idx;
  double peob, peoe;
  int N = 0;
  b = 0;
  // set_wall(wall,0);
  // set_wall(wall,plen);
  for (int i = 1; i <= plen; i++)
    if (i == plen || is_error(wall,i-1) != is_error(wall,i) || (!is_error(wall,i) && is_wall_by(OTHERS,wall,i)))
      { // set_wall(wall,i);
        e = i;
        _idx = bs_eintvl(eintvl,0,NS,b,e);
        intvl[N].b = b;
        intvl[N].e = e;
        intvl[N].cb = profile[b];
        intvl[N].ce = profile[e-1];
        intvl[N].is_rel = false;
        intvl[N].pe = (_idx != -1) ? log(eintvl[_idx].pe) : -INFINITY;
        peob = MAX(perror[b][OTHERS][DROP],perror[b][OTHERS][GAIN]);
        peoe = MAX(perror[e][OTHERS][DROP],perror[e][OTHERS][GAIN]);
        intvl[N].pe_o.b = (peob != -INFINITY) ? log(peob) : -INFINITY;
        intvl[N].pe_o.e = (peoe != -INFINITY) ? log(peoe) : -INFINITY;
        intvl[N].asgn = N_STATE;   // unclassified
        N++;
        b = e;
      }

#ifdef DEBUG_WALL
  fprintf(stderr,"Intvls:\n");
  for (int i = 0; i < N; i++)
    fprintf(stderr,"    (%d, %d) pe_S=%lf, pe_O=(%lf, %lf)\n",
                   intvl[i].b,intvl[i].e,intvl[i].pe,intvl[i].pe_o.b,intvl[i].pe_o.e);
#endif

  return N;
}

static void correct_wall_cnt(Intvl *intvl, int i, const uint16 *profile, Seq_Ctx *ctx[N_WTYPE], const int K)
{ Intvl I = intvl[i];
  pos_t first, last;
  int n_gain = 0, n_drop = 0;
  int lmax;

  last = MIN(I.b+K-1,I.e-1);
  for (pos_t i = I.b; i < last; i++)
    n_gain += MAX((int)profile[i+1]-profile[i],0);
  if (I.b+K-1 < I.e)
    { lmax = 0;
      for (int t = HP; t <= TS; t++)
        { int l = ctx[GAIN][I.b+K-1][t]*(t+1);
          if (lmax < l)
            lmax = l;
        }
      last = I.b+lmax;
      for (pos_t i = I.b; i < last; i++)
        n_gain -= MAX((int)profile[i]-profile[i+1],0);
    }
  
  first = MAX(I.e-K+1,I.b);
  for (pos_t i = first; i < I.e-1; i++)
    n_drop += MAX((int)profile[i]-profile[i+1],0);
  if (I.b < I.e-K+1)
    { lmax = 0;
      for (int t = HP; t <= TS; t++)
        { int l = ctx[DROP][I.e-K+1][t]*(t+1);
          if (lmax < l)
            lmax = l;
        }
      first = I.e-lmax;
      for (pos_t i = first; i < I.e-1; i++)
        n_drop -= MAX((int)profile[i+1]-profile[i],0);
    }

  intvl[i].ccb = MIN(I.cb+MAX(n_gain,0),MAX_KMER_CNT);
  intvl[i].cce = MIN(I.ce+MAX(n_drop,0),MAX_KMER_CNT);

  last = MIN(I.b+2*K,I.e);
  for (pos_t i = I.b; i < last; i++)
    if (intvl[i].ccb < profile[i])
      intvl[i].ccb = profile[i];
  first = MAX(I.e-2*K,I.b);
  for (pos_t i = first; i < I.e; i++)
    if (intvl[i].cce < profile[i])
      intvl[i].cce = profile[i];

#ifdef DEBUG_COR
  fprintf(stderr,"@ (%d,%d-1) %d(=%d+%d), %d(=%d+%d)\n",
                 I.b,I.e,intvl[i].ccb,I.cb,n_gain,intvl[i].cce,I.ce,n_drop);
#endif

  return;
}

int find_rel_intvl(Intvl *intvl, int N, Intvl *rintvl, cnt_t *profile, Seq_Ctx *ctx[N_WTYPE], int K)
{ int M = 0;
  double logpthres = log(PE_THRES[FINAL][SELF]);
  for (int i = 0; i < N; i++)
    { // Reliable intervals must be long and not explained by errors in others
      if (intvl[i].e-intvl[i].b < (pos_t)K)
        continue;
      if (MAX(intvl[i].cb,intvl[i].ce) >= GLOBAL_COV[REPEAT])
        continue;
      if (intvl[i].pe >= logpthres)
        continue;
      correct_wall_cnt(intvl,i,profile,ctx,K);
      if (logp_trans(intvl[i].b,intvl[i].e,intvl[i].ccb,intvl[i].cce,
                     (intvl[i].ccb+intvl[i].cce)/2) < THRES_DIFF_REL)
        continue;
      intvl[i].is_rel = true;
      rintvl[M] = intvl[i];
      M++;
    }
#ifdef DEBUG_COR
  fprintf(stderr,"Intvls (#=%d):\n",N);
  for (int i = 0; i < N; i++)
    fprintf(stderr,"    (%d, %d): c=(%d, %d), cc=(%d, %d), pe_S=%lf, pe_O=(%lf, %lf)\n",
                   intvl[i].b,intvl[i].e,intvl[i].cb,intvl[i].ce,intvl[i].ccb,intvl[i].cce,
                   intvl[i].pe,intvl[i].pe_o.b,intvl[i].pe_o.e);
  fprintf(stderr,"Rel Intvls (#=%d):\n",M);
  for (int i = 0; i < M; i++)
    fprintf(stderr,"    (%d, %d): c=(%d, %d), cc=(%d, %d), pe_S=%lf, pe_O=(%lf, %lf)\n",
                   rintvl[i].b,rintvl[i].e,rintvl[i].cb,rintvl[i].ce,rintvl[i].ccb,rintvl[i].cce,
                   rintvl[i].pe,rintvl[i].pe_o.b,rintvl[i].pe_o.e);
#endif
  return M;
}
