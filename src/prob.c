#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "ClassPro.h"
// #include <gsl/gsl_sf_bessel.h>   // TODO: use GSL instead (underflow error)
#include "bessel.h"
#include "bessel.c"

double logfact[MAX_KMER_CNT+1];

static inline void precompute_logfact()
{ for (cnt_t n = 1; n <= MAX_KMER_CNT; n++)
    logfact[n] = logfact[n-1]+log(n);

  return;
}

#ifdef DEBUG
static inline int _check_cnt(cnt_t n)
{
  if (n > MAX_KMER_CNT)
    { fprintf(stderr,"K-mer count (%d) > MAX_KMER_CNT (%d) (due to D/R ratio?)\n",n,MAX_KMER_CNT);
      // exit(1);
      return MAX_KMER_CNT;
    }
  return n;
}
#endif

static inline double logp_poisson(cnt_t k, int lambda)
{ 
#ifdef DEBUG
  k = _check_cnt(k);
#endif
  return k * log((double)lambda) - lambda - logfact[k];
}

static inline double logp_skellam(int k, double lambda)
{ return -2. * lambda + log(bessi(abs(k),2.*lambda));   // TODO: precompute bessi?
  // return -2. * lambda + log(gsl_sf_bessel_In(abs(k),2.*lambda));
}

#ifdef DEBUG
static inline void _check_cnt_binom(cnt_t *k, cnt_t *n)
{ 
  *k = _check_cnt(*k);
  *n = _check_cnt(*n);
  if (*k > *n)
    { fprintf(stderr,"k (%d) > n (%d) in Binom\n",*k,*n);
      exit(1);
    }
  return;
}
#endif

static inline double logp_binom(cnt_t k, cnt_t n, double p)
{ 
#ifdef DEBUG
  _check_cnt_binom(&k,&n);
#endif
  return logfact[n] - logfact[k] - logfact[n-k] + k * log(p) + (n-k) * log(1-p);
}

static inline double logp_binom_pre(cnt_t k, cnt_t n, double lpe, double l1mpe)
{ 
#ifdef DEBUG
  _check_cnt_binom(&k,&n);
#endif
  return logfact[n] - logfact[k] - logfact[n-k] + k * lpe + (n-k) * l1mpe;
}

// TODO: chi-square is faster when k,n are large?
static inline double binom_test_g(cnt_t k, cnt_t n, double pe, bool exact)
{ 
#ifdef DEBUG
  _check_cnt_binom(&k,&n);
#endif

  const double lpe      = log(pe);   // TODO: precompute in error models?
  const double l1mpe    = log(1-pe);
  const double mean     = n * pe;
  const bool   decrease = ((double)k >= mean) ? true : false;
  double p, p_first, p_curr;

#ifdef DEBUG_BINOM
  fprintf(stderr,"BinomTest(%d; %d, %lf) [%s]\n",k,n,pe,exact ? "exact" : "approx");
  fprintf(stderr,"    k (%d) %s mean (%lf)\n",k,decrease ? ">=" : "<",mean);
#endif

  if (decrease)
    { p = p_first = exp(logp_binom_pre(k,n,lpe,l1mpe));
      for (cnt_t x = k+1; x <= n; x++)
        { p += p_curr = exp(logp_binom_pre(x,n,lpe,l1mpe));
          if (!exact && 10 * p_curr < p_first)
            break;
        }
    }
  else
    { p = p_first = (k == 0) ? 0. : exp(logp_binom_pre(k-1,n,lpe,l1mpe));
      for (int x = k-2; x >= 0; x--)
        { p += p_curr = exp(logp_binom_pre(x,n,lpe,l1mpe));
          if (!exact && 10 * p_curr < p_first)
            break;
        }
      p = 1-p;
    }

  return p;
}
