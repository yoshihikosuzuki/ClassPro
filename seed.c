#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "nthash.h"
#include "kdq.h"
#include "ClassPro.h"

#define DEBUG_HASH
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

static void kmer_hash(const char *seq, int *hash, int plen, int K)
{ uint64_t hVal, fhVal = 0, rhVal = 0;   // canonical, forward, and reverse-strand hash values
  hVal = NTC64_b(seq-K+1,K,&fhVal,&rhVal);   // initial hash value
  hash[0] = (int)(hVal % MOD);   // TODO: no need of taking MOD? (a little space efficient, though)
#ifdef DEBUG_HASH
  fprintf(stderr,"Canonical hash values:\n");
  fprintf(stderr,"  i = %2d: kmer = %.*s, hash = %d (h = %lu, f = %lu, r = %lu)\n",
                 0,K,seq-K+1,hash[0],hVal,fhVal,rhVal);
#endif

  for (int i = 1; i < plen; i++)
    { hVal = NTC64_c(seq[i-K],seq[i],K,&fhVal,&rhVal); // consecutive hash values
      hash[i] = (int)(hVal % MOD);
#ifdef DEBUG_HASH
      fprintf(stderr,"  i = %2d: kmer = %.*s, hash = %d (h = %lu, f = %lu, r = %lu)\n",
                     i,K,seq-K+1+i,hash[i],hVal,fhVal,rhVal);
#endif
    }
  return;
}

static void _find_seeds(const char *seq, const uint16 *profile, const char *class, const int *hash, const int plen, const int K, int *sasgn, const char C)
{ 
#if defined(DEBUG_SEED) || defined(INFO_SEED)
  fprintf(stderr,"\n");
#endif

  kdq_t(hmer_t) *Q = kdq_init(hmer_t);
  hmer_t e, f;
  hmer_t *p;
  
  // Uniformly sparse count maximizers using sliding window maximum alrogithm
  // sasgn:   -10 = not candidate,   -2 = seed fixed elsewhere,   -1 = fixed seed,   0 = candidate,   1 = seed
  for (int i = 0; i < plen; i++)
    if (sasgn[i] != -2)
      sasgn[i] = (class[i] == C) ? 0 : -10;
  bool converged = false;
  int iter = 1;
  while (!converged)
    { for (int i = 0; i < plen; i++)   // TODO: Currently using deque. How about binomial heap?
        { if (sasgn[i] >= 0)
            { // Add a new element `e` into the deque
              e.pos = i;
              e.key = profile[i];
              while (kdq_size(Q) > 0) {
                f = kdq_at(Q, kdq_size(Q) - 1);   // TODO: return pointer instead of object?
                if (f.key < e.key) {
                  kdq_size(Q)--;
                } else {
                  break;
                }
              }
              kdq_push(hmer_t, Q, e);
            }
          if (kdq_size(Q) == 0) continue;
          // Remove out-of-range elements
          if (kdq_at(Q, 0).pos <= i - WSIZE) {
            while (kdq_size(Q) > 0 && kdq_at(Q, 0).pos <= i - WSIZE) {
              p = kdq_shift(hmer_t, Q);
            }
          }
          if (kdq_size(Q) == 0) continue;
          // Now the leftmost element in the deque is the maximizer in the current window.
          // Increament `sasgn` so that `sasgn[i]` := # of windows in which k-mer at `i` is the maximizer.
          e = kdq_at(Q, 0);
          sasgn[e.pos]++;
          // Increment for each k-mer with tie count   // TODO: Worst case O(L^2). Make it efficient
          for (int j = 1; j < (int)kdq_size(Q); j++) {
            f = kdq_at(Q, j);
            if (e.key != f.key) break;
            sasgn[f.pos]++;
          }
#ifdef DEBUG_SEED
          // fprintf(stderr, "@ %5d (%c): |Q| = %ld, Q =", i, class[i], kdq_size(Q));
          // for (int j = 0; j < (int)kdq_size(Q); j++) {
          //   fprintf(stderr, "  %d @ %d (%d)", kdq_at(Q, j).key, kdq_at(Q, j).pos, sasgn[kdq_at(Q, j).pos]);
          // }
          // fprintf(stderr, "\n");
#endif
        }   // for (int i = 0; i < plen - WSIZE; i++)
#ifdef DEBUG_SEED
      for (int i = 0; i < plen; i++)
        if (sasgn[i] >= 0)
          { if (class[i] != C)
              { fprintf(stderr, "sasgn[%d] >= 0 where class(%c) != %c\n", i, class[i], C);
                exit(1);
              }
            fprintf(stderr,"%c-maximizer(%d) @ %5d: kmer = %.*s, count = %d, #window = %d\n",C,iter,i,K,seq+i-K+1,profile[i],sasgn[i]);
          }
#endif
      // Find USM
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == WSIZE) sasgn[i] = -1;
      // Filter out interval (i-|w|,i+|W|) for each position i where the k-mer at i is USM
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == -1)
          { for (int j = MAX(0,i-WSIZE+1); j < MIN(plen,i+WSIZE); j++)   // FIXME: Make it efficient
              if (sasgn[j] >= 0)
                sasgn[j] = -10;
          }
      converged = true;
      for (int i = 0; i < plen; i++)
        if (sasgn[i] >= 0)
          { converged = false;
            sasgn[i] = 0;
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

  exit(0);

  // Uniformly sparse modimizers among the maximizers
  iter = 1; converged = false;
  while (!converged)
    { for (int i = 0; i < plen-WSIZE; i++)
        { int hmin = MOD;
          for (int j = 0; j < WSIZE; j++)
            if (sasgn[i+j] >= 0) hmin = MIN(hash[i+j],hmin);
          for (int j = 0; j < WSIZE; j++)
            if (sasgn[i+j] == 1 && hash[i+j] > hmin)
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
void find_seeds(const char *_seq, const uint16 *profile, const char *class, const int plen, const int K, int *sasgn, int *hash)
{ const char *seq = _seq+K-1;   // kmer seq @ i on profile = seq[i-K+1]..seq[i]
  // Compute canonical hash for every k-mer
  kmer_hash(seq,hash,plen,K);

  // Find seeds from first H-mers and then D-mers
  _find_seeds(seq,profile,class,hash,plen,K,sasgn,'H');
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == 1) sasgn[i] = -2;
  _find_seeds(seq,profile,class,hash,plen,K,sasgn,'D');
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
