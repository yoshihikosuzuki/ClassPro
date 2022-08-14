#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "nthash.h"
#include "ClassPro.h"

#undef DEBUG_HASH
#undef DEBUG_SEED
#undef INFO_SEED

#define WSIZE 1000
#define MOD 1009

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

static inline void count_maximizer(const char *seq, const uint16 *profile, const char *class, const int plen, const int K, int *sasgn, const char C, kdq_t(hmer_t) *Q)
{ //kdq_t(hmer_t) *Q = kdq_init(hmer_t);   // TODO: reuse deque throughout CP?
  hmer_t e, f;
  hmer_t *p;

#ifdef DEBUG_SEED
  if (kdq_size(Q) != 0)
    { fprintf(stderr,"deque not empty!\n");
      exit(1);
    }
#endif

  // Uniformly sparse count maximizers using sliding window maximum algorithm
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
        }   // for (int i = 0; i < plen; i++)

#ifdef DEBUG_SEED
      for (int i = 0; i < plen; i++)
        if (sasgn[i] >= 0)
          { if (class[i] != C)
              { fprintf(stderr, "sasgn[%d] >= 0 where class(%c) != %c\n", i, class[i], C);
                exit(1);
              }
            bool is_usm = false;
            if (   (i < plen-WSIZE+1 && sasgn[i] == WSIZE)
                || (plen-WSIZE+1 <= i && sasgn[i] == plen-i) )
                { is_usm = true; }
            fprintf(stderr,"%c-maximizer(%d) %c @ %5d: kmer = %.*s, count = %2d, #window = %4d\n",
                           C,iter,(is_usm) ? '*' : ' ',i,K,seq+i-K+1,profile[i],sasgn[i]);
          }
#endif

      // Find USM (while being careful for both `WSIZE`-bp ends of a read)
      for (int i = 0; i < plen-WSIZE+1; i++)
        if (sasgn[i] == WSIZE) sasgn[i] = -1;
      for (int i = plen-WSIZE+1; i < plen; i++)
        if (sasgn[i] == plen-i) sasgn[i] = -1;

      // Filter out interval (i-|w|,i+|W|) for each position i where the k-mer at i is USM
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == -1)
          { for (int j = MAX(0,i-WSIZE+1); j < MIN(plen,i+WSIZE); j++)   // FIXME: Make it efficient
              if (sasgn[j] >= 0)
                sasgn[j] = -10;
          }

      // Converged?
      converged = true;
      for (int i = 0; i < plen; i++)
        if (sasgn[i] >= 0)
          { converged = false;
            sasgn[i] = 0;
          }
      // Empty the queue   // TODO: which of the two below is better?
      // while (kdq_pop(hmer_t,Q)) {}
      kdq_size(Q) = 0;
      iter++;

#ifdef DEBUG_SEED
      // if (iter > 10) break;
#endif
    }   // while (!converged)

  // Change flag to "candidate" for the next round, i.e. minimizer
#ifdef INFO_SEED
  int c = 0;
#endif
  for (int i = 0; i < plen; i++)
    {
#ifdef DEBUG_SEED
      if (sasgn[i] >= 0)
        { fprintf(stderr,"sasgn[%d] = %d >= 0 after convergence!\n",i,sasgn[i]);
          exit(1);
        }
#endif
      if (sasgn[i] == -1)
        { sasgn[i] = 0;
#ifdef DEBUG_SEED
          fprintf(stderr,"%c-maximizer @ %5d: kmer = %.*s, count = %2d\n",C,i,K,seq+i-K+1,profile[i]);
#endif
#ifdef INFO_SEED
          c++;
#endif
        }
    }
#ifdef INFO_SEED
  fprintf(stderr,"%d %c-maximizers\n",c,C);
#endif

  // kdq_destroy(hmer_t,Q);

  return;
}

static inline void hash_minimizer(const char *seq, const uint16 *profile, const char *class, const int *hash, const int plen, const int K, int *sasgn, const char C, kdq_t(hmer_t) *Q)
{ // kdq_t(hmer_t) *Q = kdq_init(hmer_t);
  hmer_t e, f;
  hmer_t *p;

#ifdef DEBUG_SEED
  if (kdq_size(Q) != 0)
    { fprintf(stderr,"deque not empty!\n");
      exit(1);
    }
#endif

  // Uniformly sparse modimizers among the maximizers
  bool converged = false;
  int iter = 1;
  while (!converged)
    { for (int i = 0; i < plen; i++)
        { if (sasgn[i] >= 0)
            { e.pos = i;
              e.key = hash[i];
              while (kdq_size(Q) > 0) {
                f = kdq_at(Q, kdq_size(Q)-1);
                if (f.key > e.key) {
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
          // Now the leftmost element in the deque is the minimizer in the current window.
          // Increament `sasgn` so that `sasgn[i]` := # of windows in which k-mer at `i` is the maximizer.
          e = kdq_at(Q, 0);
          sasgn[e.pos]++;
          // Increment for each k-mer with tie hash   // TODO: Worst case O(L^2). Make it efficient
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
        }   // for (int i = 0; i < plen; i++)

#ifdef DEBUG_SEED
      for (int i = 0; i < plen; i++)
        if (sasgn[i] >= 0)
          { if (class[i] != C)
              { fprintf(stderr, "sasgn[%d] >= 0 where class(%c) != %c\n", i, class[i], C);
                exit(1);
              }
            bool is_usm = false;
            if (   (i < plen-WSIZE+1 && sasgn[i] == WSIZE)
                || (plen-WSIZE+1 <= i && sasgn[i] == plen-i) )
                { is_usm = true; }
            fprintf(stderr,"%c-minimizer(%d) %c @ %5d: kmer = %.*s, count = %2d, hash = %5d, #window = %4d\n",
                           C,iter,(is_usm) ? '*' : ' ',i,K,seq+i-K+1,profile[i],hash[i],sasgn[i]);
          }
#endif

      // Find USM (while being careful for both `WSIZE`-bp ends of a read)
      for (int i = 0; i < plen-WSIZE+1; i++)
        if (sasgn[i] == WSIZE) sasgn[i] = -1;
      for (int i = plen-WSIZE+1; i < plen; i++)
        if (sasgn[i] == plen-i) sasgn[i] = -1;

      // Filter out interval (i-|w|,i+|W|) for each position i where the k-mer at i is USM
      for (int i = 0; i < plen; i++)
        if (sasgn[i] == -1)
          { for (int j = MAX(0,i-WSIZE+1); j < MIN(plen,i+WSIZE); j++)   // FIXME: Make it efficient
              if (sasgn[j] >= 0)
                sasgn[j] = -10;
          }

      // Converged?
      converged = true;
      for (int i = 0; i < plen; i++)
        if (sasgn[i] >= 0)
          { converged = false;
            sasgn[i] = 0;
          }
      // Empty the queue   // TODO: which of the two below is better?
      // while (kdq_pop(hmer_t,Q)) {}
      kdq_size(Q) = 0;
      iter++;

#ifdef DEBUG_SEED
      // if (iter > 10) break;
#endif
    }   // while (!converged)

  // Change flag to "fixed seed"
#ifdef INFO_SEED
  int c = 0;
#endif
  for (int i = 0; i < plen; i++)
    {
#ifdef DEBUG_SEED
      if (sasgn[i] >= 0)
        { fprintf(stderr,"sasgn[%d] = %d >= 0 after convergence!\n",i,sasgn[i]);
          exit(1);
        }
#endif
      if (sasgn[i] == -1)
        { sasgn[i] = -2;
#ifdef DEBUG_SEED
          fprintf(stderr,"%c-minimizer @ %5d: kmer = %.*s, count = %2d, hash = %5d\n",C,i,K,seq+i-K+1,profile[i],hash[i]);
#endif
#ifdef INFO_SEED
          c++;
#endif
        }
    }
#ifdef INFO_SEED
  fprintf(stderr,"%d %c-minimizers\n",c,C);
#endif

  // kdq_destroy(hmer_t,Q);

  return;
}

static void _find_seeds(const char *seq, const uint16 *profile, const char *class, const int *hash, const int plen, const int K, int *sasgn, const char C, kdq_t(hmer_t) *Q)
{ // Meaning of the value of `sasgn` (for each position (= k-mer)):   // TODO: change to bit operation with macros?
  //   -10 = not candidate
  //    -2 = seed fixed elsewhere
  //    -1 = fixed seed
  //     0 = candidate
  //     1 = seed
  for (int i = 0; i < plen; i++)
    if (sasgn[i] != -2)
      sasgn[i] = (class[i] == C) ? 0 : -10;

  count_maximizer(seq,profile,class,plen,K,sasgn,C,Q);
  hash_minimizer(seq,profile,class,hash,plen,K,sasgn,C,Q);

  return;
}

// NOTE: len(profile) == len(class) == plen, len(_seq) = plen + K - 1
// NOTE: class[i] in {'E', 'H', 'D', 'R'}
void find_seeds(const char *_seq, const uint16 *profile, const char *class, const int plen, const int K, int *sasgn, int *hash, kdq_t(hmer_t) *Q)
{
#if defined(DEBUG_SEED) || defined(INFO_SEED)
  fprintf(stderr,"\n");
#endif

  // printf("c");
  // for (int i = -(K-1); i < plen; i++)
  //   printf("%c",class[i]);
  // printf("\n");

  const char *seq = _seq+K-1;   // kmer seq @ i on profile = seq[i-K+1]..seq[i]
  // Compute canonical hash for every k-mer
  kmer_hash(seq,hash,plen,K);

  // Find seeds from first H-mers and then D-mers
  _find_seeds(seq,profile,class,hash,plen,K,sasgn,'H',Q);
  _find_seeds(seq,profile,class,hash,plen,K,sasgn,'D',Q);   // TODO: density adjustment based on the H-MMs

  // change flag value
  for (int i = 0; i < plen; i++)
    { // sasgn[i] = (sasgn[i] == -2) ? class[i] : 'E';
      if (sasgn[i] == -2)   // H/D-seed
        sasgn[i] = class[i];
      else if (class[i] == 'E')
        sasgn[i] = 'E';
      else
        sasgn[i] = 'R';
    }

  // printf("s");
  // for (int i = 0; i < K-1; i++)
  //   printf(" ");
  // for (int i = 0; i < plen; i++)
  //   printf("%c",sasgn[i]);
  // printf("\n");

#ifdef DEBUG_SEED
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == 'H' || sasgn[i] == 'D')
      fprintf(stderr,"seed(%c) @ %5d: kmer = %.*s, count = %d\n",class[i],i,K,seq+i-K+1,profile[i]);
#endif

#ifdef DEBUG_ITER
  { int c = 0, e = 0;
    for (int i = 0; i < plen; i++)
      { if (sasgn[i] == 'H' || sasgn[i] == 'D') c++;
        else if (sasgn[i] == 'E') e++;
      }
    fprintf(stderr,", %3d seeds, %3d errors",c,e);
  }
#endif
  return;
}
