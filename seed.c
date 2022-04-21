#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "kdq.h"
#include "ClassPro.h"

#define DEBUG_SEED
#define INFO_SEED

#define WSIZE 1000
#define MOD 1009

static const char ntoi[128] = { 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 1, 0, 0, 0, 2,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 3, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0  };
static const char ntor[128] = { 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 3, 0, 2, 0, 0, 0, 1,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0  };

typedef struct {
  int pos;
  int key;
} hmer_t;

KDQ_INIT(hmer_t);

static inline bool is_canonical(const char *str, int K)
{ for (int i = 0; i < K/2+K%2; i++)
    if (ntoi[(int)str[i]] != (ntor[(int)str[K-1-i]]))
      return (ntoi[(int)str[i]] < (ntor[(int)str[K-1-i]]));
#ifdef DEBUG
  if (K%2 == 1)
    { fprintf(stderr,"K is odd but palindrome\n");
      exit(1);
    }
#endif
  return false;   // palindrome
}

static inline char complement(char c)
{ switch (c)
    { case 'A':
        return 'T';
      case 'C':
        return 'G';
      case 'G':
        return 'C';
      case 'T':
        return 'A';
      default:
        fprintf(stderr,"Not a nucleotide\n");
        exit(1);
    }
}

static inline int kmer_hash(const char *str, int K)
{ unsigned long hash = 5381;
  int c;
  bool is_can = is_canonical(str,K);
  for (int i = 0; i < K; i++)
    { c = (is_can) ? str[i] : complement(str[K-1-i]);
      hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }
  return (int)(hash % MOD);
}

static void _find_seeds(const char *seq, const uint16 *profile, const char *class, const int plen, const int K, int *sasgn, const char C)
{ 
#if defined(DEBUG_SEED) || defined(INFO_SEED)
  fprintf(stderr,"\n");
#endif

  kdq_t(hmer_t) *Q = kdq_init(hmer_t);
  hmer_t e;
  hmer_t *p;
  
  // Uniformly sparse count maximizers
  // sasgn = -10: not candidate, -2: seed fixed elsewhere, -1: fixed seed, 0: candidate, 1: seed
  for (int i = 0; i < plen; i++)
    if (sasgn[i] != -2)
      sasgn[i] = (class[i] == C) ? 0 : -10;
  bool converged = false;
  int iter = 1;
  while (!converged)
    { for (int i = 0; i < plen - WSIZE; i++)   // TODO: deque or binomial heap
        { if (sasgn[i] >= 0)
            { e.pos = i;
              e.key = profile[i];
              while (kdq_size(Q) > 0) {
                hmer_t f = kdq_at(Q, kdq_size(Q) - 1);   // TODO: return pointer instead of object?
                if (f.key < e.key) {
                  kdq_size(Q)--;
                } else {
                  break;
                }
              }
              kdq_push(hmer_t, Q, e);
            }
          if (kdq_size(Q) == 0) continue;
          if (kdq_at(Q, 0).pos <= i - WSIZE) {
            while (kdq_size(Q) > 0 && kdq_at(Q, 0).pos <= i - WSIZE) {
              p = kdq_shift(hmer_t, Q);
            }
          }
          if (kdq_size(Q) == 0) continue;
          e = kdq_at(Q, 0);   // TODO: tie?
          sasgn[e.pos]++;
          // int cmax = -1;
          // for (int j = 0; j < WSIZE; j++)
          //   if (sasgn[i+j] >= 0) cmax = MAX(profile[i+j],cmax);
          // for (int j = 0; j < WSIZE; j++)
          //   if (sasgn[i+j] == 1 && profile[i+j] < cmax)
          //     sasgn[i+j] = 0;
#ifdef DEBUG_SEED
          // fprintf(stderr, "@ %5d: |Q| = %ld, Q =", i, kdq_size(Q));
          // for (int j = 0; j < (int)kdq_size(Q); j++) {
          //   fprintf(stderr, "  %d @ %d", kdq_at(Q, j).key, kdq_at(Q, j).pos);
          // }
          // fprintf(stderr, "\n");
#endif
        }
#ifdef DEBUG_SEED
      // for (int i = 0; i < plen; i++)
      //   if (sasgn[i] >= 0) fprintf(stderr,"%c-maximizer(%d) @ %5d: kmer = %.*s, count = %d, sasgn = %d\n",C,iter,i,K,seq+i-K+1,profile[i], sasgn[i]);
      // exit(1);
#endif
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == 1)
          { sasgn[i] = -1;
            for (int j = MAX(0,i-WSIZE+1); j < MIN(plen,i+WSIZE); j++)
              if (sasgn[j] == 0)
                sasgn[j] = -10;
          }
      converged = true;
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == 0)
          { converged = false;
            sasgn[i] = 1;
          }
      iter++;
    }
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == -1) sasgn[i] = 1;
#ifdef DEBUG_SEED
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == 1) fprintf(stderr,"%c-maximizer @ %5d: kmer = %.*s, count = %d\n",C,i,K,seq+i-K+1,profile[i]);
#endif

#ifdef INFO_SEED
  { int c = 0;
    for (int i = 0; i < plen; i++)
      if (sasgn[i] == 1) c++;
    fprintf(stderr,"%d %c-maximizers\n",c,C);
  }
#endif

  // Uniformly sparse modimizers among the maximizers
  iter = 1; converged = false;
  while (!converged)
    { for (int i = 0; i < plen-WSIZE; i++)
        { int hmin = MOD;
          for (int j = 0; j < WSIZE; j++)
            if (sasgn[i+j] >= 0) hmin = MIN(kmer_hash(seq+i+j-K+1,K),hmin);
          for (int j = 0; j < WSIZE; j++)
            if (sasgn[i+j] == 1 && kmer_hash(seq+i+j-K+1,K) > hmin)
              sasgn[i+j] = 0;
        }
#ifdef DEBUG_SEED
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == 1) fprintf(stderr,"%c-minimizer(%d) @ %5d: kmer = %.*s, count = %d\n",C,iter,i,K,seq+i-K+1,profile[i]);
#endif
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == 1)
          { sasgn[i] = -1;
            for (int j = MAX(0,i-WSIZE+1); j < MIN(plen,i+WSIZE); j++)
              if (sasgn[j] == 0)
                sasgn[j] = -10;
          }
      converged = true;
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == 0)
          { converged = false;
            sasgn[i] = 1;
          }
      iter++;
    }
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == -1) sasgn[i] = 1;
#ifdef DEBUG_SEED
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == 1) fprintf(stderr,"%c-minimizer @ %5d: kmer = %.*s, count = %d\n",C,i,K,seq+i-K+1,profile[i]);
#endif

#ifdef INFO_SEED
  { int c = 0;
    for (int i = 0; i < plen; i++)
      if (sasgn[i] == 1) c++;
    fprintf(stderr,"%d %c-minimizers\n",c,C);
  }
#endif

  return;
}

// NOTE: len(profile) == len(class) == plen, len(_seq) = plen + K - 1
// NOTE: class[i] in {'E', 'H', 'D', 'R'}
void find_seeds(const char *_seq, const uint16 *profile, const char *class, const int plen, const int K, int *sasgn)
{ const char *seq = _seq+K-1;   // kmer seq @ i on profile = seq[i-K+1]..seq[i]
  _find_seeds(seq,profile,class,plen,K,sasgn,'H');
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == 1) sasgn[i] = -2;
  _find_seeds(seq,profile,class,plen,K,sasgn,'D');
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == 1) sasgn[i] = -2;

  // change flag value
  for (int i = 0; i < plen; i++)
    sasgn[i] = (sasgn[i] == -2) ? class[i] : 'E';

#ifdef DEBUG_SEED
  for (int i = 0; i < plen; i++)
    if (sasgn[i] != 'E') fprintf(stderr,"seed(%c) @ %5d: kmer = %.*s, count = %d\n",class[i],i,K,seq+i-K+1,profile[i]);
#endif

#ifdef INFO_SEED
  { int c = 0;
    for (int i = 0; i < plen; i++)
      if (sasgn[i] != 'E') c++;
    fprintf(stderr,"%d seeds\n",c);
  }
#endif
  return; 
}
