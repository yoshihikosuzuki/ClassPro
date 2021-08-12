#ifndef _PROB_H
#define _PROB_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "ClassPro.h"
#include "bessel.h"

extern double logfact[32768];

void precompute_logfact();

static inline double logp_poisson(int k, int lambda)
{ k = MIN(k,32767);   // TODO: * dr_ratio can make imaginary count > 32767 (but should not be such a large count, why?)
  return k * log((double)lambda) - lambda - logfact[k];
}

static inline double logp_skellam(int k, double lambda)
{ return -2. * lambda + log(bessi(abs(k),2.*lambda));   // TODO: precompute bessi?
}

static inline double logp_binom(int k, int n, double p)
{ k = MIN(k,32767);
  n = MIN(n,32767);
  return logfact[n] - logfact[k] - logfact[n-k] + k * log(p) + (n-k) * log(1-p);
}

static inline double logp_binom_pre(int k, int n, double lpe, double l1mpe)
{ k = MIN(k,32767);
  n = MIN(n,32767);
  return logfact[n] - logfact[k] - logfact[n-k] + k * lpe + (n-k) * l1mpe;
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

#endif // _PROB_H