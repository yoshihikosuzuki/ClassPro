#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "nthash.h"
#include "ClassPro.h"

#undef DEBUG_HASH
#undef DEBUG_COMPRESS
#undef DEBUG_QUEUE
#undef DEBUG_SEED
#undef DEBUG_SEGMENT
#undef DEBUG_SELECT
#undef TIME_SEED
#undef INFO_SEED

#define WSIZE 1000
#define MOD 10009

static void kmer_hash(const char *seq,
                      int        *hash,
                      int         plen,
                      int         K)
{ uint64_t hVal, fhVal, rhVal;   // canonical, forward, and reverse-strand hash values

  fhVal = rhVal = 0;
  hVal = NTC64_b(seq-K+1,K,&fhVal,&rhVal);   // initial hash value
  hash[0] = (int)(hVal % MOD);   // TODO: no need of taking MOD? (a little space efficient, though)

#ifdef DEBUG_HASH
  fprintf(stderr,"Canonical hash values:\n");
  fprintf(stderr,"  i = %2d: kmer = %.*s, hash = %d (h = %lu, f = %lu, r = %lu)\n",
                 0,K,seq-K+1,hash[0],hVal,fhVal,rhVal);
#endif

  for (int i = 1; i < plen; i++)
    { hVal = NTC64_c(seq[i-K],seq[i],K,&fhVal,&rhVal);   // consecutive hash values
      hash[i] = (int)(hVal % MOD);

#ifdef DEBUG_HASH
      fprintf(stderr,"  i = %2d: kmer = %.*s, hash = %d (h = %lu, f = %lu, r = %lu)\n",
                     i,K,seq-K+1+i,hash[i],hVal,fhVal,rhVal);
#endif
    }

  return;
}

static int compress_profile(const uint16 *profile,
                            const char   *class,
                            seg_t        *cprofile,
                            const int     plen,
                            const char    C)
{ int N = 0;
  int b = 0, e = 1;
  bool prev_valid = (class[0] == C) ? true : false;

  while (e < plen)
    { if (!prev_valid)
        { while (e < plen && class[e] != C)
            e++;
          cprofile[N].b = b;
          cprofile[N].e = e;
          cprofile[N].cnt = -1;
          cprofile[N].nw = -10;
          cprofile[N].is_seed = false;
          N++;
          b = e;
          e++;
          prev_valid = true;
        }
      else
        { while (e < plen && profile[e] == profile[e-1])
            e++;
          cprofile[N].b = b;
          cprofile[N].e = e;
          cprofile[N].cnt = profile[e-1];
          cprofile[N].nw = 0;
          cprofile[N].is_seed = false;
          N++;
          b = e;
          e++;
          prev_valid = (class[b] == C) ? true : false;
        }
    }

#ifdef DEBUG_COMPRESS
  fprintf(stderr,"# of segments = %d\n",N);
  for (int i = 0; i < N; i++)
    fprintf(stderr,"%d: (%d, %d) cnt = %d\n",i,cprofile[i].b,cprofile[i].e,cprofile[i].cnt);
#endif

  return N;
}

static int compare_cprofile(const void *a, const void *b)
{ return ((seg_t*)b)->nw - ((seg_t*)a)->nw;
}

static int compare_intvl(const void *a, const void *b)
{ return ((intvl_t*)a)->b - ((intvl_t*)b)->b;
}

static inline bool does_ovlp(int ab, int ae, int bb, int be)
{ return (MAX(ab, bb) <= MIN(ae-1, be-1));
}

static int bs_mintvl(intvl_t *intvl, int l, int r, int b, int e) {
  if (l > r)
     return -1;
  int m = (l+r)/2;
  if (does_ovlp(intvl[m].b,intvl[m].e,b,e))
    return m;
  else if (intvl[m].b < b)
    return bs_mintvl(intvl,m+1,r,b,e);
  else
    return bs_mintvl(intvl,l,m-1,b,e);
}

/* NOTE: Assuming interval X is always contained if there is an interval Y
         in `mintvl` that overlaps with X.
*/
static inline bool is_contained(intvl_t *mintvl, int M, int b, int e)
{ int idx = bs_mintvl(mintvl,0,M,b,e);
  if (idx == -1)
    {
#ifdef DEBUG_SELECT
      fprintf(stderr,"Not contained\n");
#endif
      return false;
    }
  else
    {
#ifdef DEBUG_SELECT
      fprintf(stderr,"Contained in mintvl[%d] = (%d, %d)\n",idx,mintvl[idx].b,mintvl[idx].e);
#endif
      return (mintvl[idx].b <= b && e <= mintvl[idx].e);
    }
}

/* Add [b..e) into `mintvl` and merge overlapping intervals
*/
static inline int add_intvl(intvl_t *mintvl, int M, int b, int e)
{ int idx = bs_mintvl(mintvl,0,M,b,e);
  if (idx == -1)
    { M++;
      mintvl[M].b = b;
      mintvl[M].e = e;
      qsort(mintvl,M,sizeof(intvl_t),compare_intvl);
      return M;
    }
  int l = idx - 1;
  while (l >= 0 && does_ovlp(mintvl[l].b,mintvl[l].e,b,e))
    l--;
  l++;
  int r = idx + 1;
  while (r < M && does_ovlp(mintvl[r].b,mintvl[r].e,b,e))
    r++;
  r--;
  mintvl[l].b = MIN(mintvl[l].b,b);
  mintvl[l].e = MAX(mintvl[r].e,e);
  if (l == r)
    return M;
  int d = r-l;
  M -= d;
  for (int i = l+1; i < M; i++)
    { mintvl[i].b = mintvl[i+d].b;
      mintvl[i].e = mintvl[i+d].e;
    }
  return M;
}

/* Meaning of the value of `sasgn` (for each position (= k-mer)):   // TODO: change to bit operation with macros?
     -10 = not candidate
      -2 = seed fixed elsewhere
      -1 = fixed seed
       0 = candidate
       1 = seed
*/
static void _find_seeds(kdq_t(hmer_t) *Q,
                        const uint16  *profile,
                        const char    *seq,
                        const char    *class,
                        const int     *hash,
                        int           *sasgn,
                        seg_t         *cprofile,
                        intvl_t       *mintvl,
                        const int      plen,
                        const int      K,
                        const char     C)
{ int N, M;
  hmer_t now, first, elem;
  seg_t seg;

  N = compress_profile(profile,class,cprofile,plen,C);

#ifdef DEBUG_SEED
  if (kdq_size(Q) != 0)
    { fprintf(stderr,"deque not empty!\n");
      exit(1);
    }
#endif

#ifdef TIME_SEED
  fprintf(stderr,"cm compress\n");
  timeTo(stderr,false);
#endif

  // Uniformly sparse count maximizers using sliding window maximum algorithm
  // and tie-count-compressed profile
  bool last_oor = false;
  int  last_oor_pos = 0;   // Position where the last Out-Of-Range happened
  for (int i = 0; i < N; i++)
    { seg = cprofile[i];
      if (seg.cnt >= 0)
        { // Add a new element `e` into the deque
          now.seg_id = i;
          now.b = seg.b;
          now.e = seg.e;
          now.cnt = seg.cnt;

          if (kdq_size(Q) > 0)
            { first = kdq_at(Q, 0);
              if (first.cnt < now.cnt)   // all elements are wiped out from Q
                { last_oor = false;
                  for (int j = 0; j < (int)kdq_size(Q); j++)
                    { elem = kdq_at(Q, j);
                      if (first.cnt == elem.cnt)
                        { cprofile[elem.seg_id].nw = MIN(now.b - elem.b,WSIZE);
#ifdef DEBUG_QUEUE
                          fprintf(stderr,"seg[%d] (cnt=%d,b=%d) nw <- %d (Found larger count)\n",
                                         elem.seg_id,elem.cnt,elem.b,cprofile[elem.seg_id].nw);
#endif
                        }
                      else
                        break;
                    }
                  kdq_size(Q) = 0;
                }
            }
          while (kdq_size(Q) > 0)
            { elem = kdq_at(Q, kdq_size(Q) - 1);
              if (elem.cnt < now.cnt)
                { kdq_size(Q)--;
                }
              else
                break;
            }
          kdq_push(hmer_t, Q, now);
        }

      if (kdq_size(Q) == 0) continue;

      // Remove out-of-range elements
      while (kdq_size(Q) > 0 && kdq_at(Q, 0).b <= seg.b - WSIZE)
        { first = kdq_at(Q, 0);
          if (last_oor)
            { cprofile[first.seg_id].nw = MIN(first.b-last_oor_pos+1,WSIZE);
#ifdef DEBUG_QUEUE
          fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Out of range & last OoR @ %d)\n",
                          first.seg_id,first.cnt,first.b,cprofile[first.seg_id].nw,last_oor_pos);
#endif
            }
          else
            { cprofile[first.seg_id].nw = WSIZE;
#ifdef DEBUG_QUEUE
          fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Out of range)\n",
                          first.seg_id,first.cnt,first.b,cprofile[first.seg_id].nw);
#endif
            }
          if (kdq_size(Q) > 1 && first.cnt > kdq_at(Q, 1).cnt)
            { last_oor_pos = first.e;
            }
          kdq_shift(hmer_t, Q);
          last_oor = true;
        }

      if (kdq_size(Q) == 0) continue;

#ifdef DEBUG_QUEUE
      fprintf(stderr,"@ %5d (%c): |Q| = %ld, Q =",i,class[i],kdq_size(Q));
      for (int j = 0; j < (int)kdq_size(Q); j++) {
        fprintf(stderr,"  %d @ seg[%d] (b=%d)",
                       kdq_at(Q, j).cnt,kdq_at(Q, j).seg_id,kdq_at(Q, j).b);
      }
      fprintf(stderr,"\n");
#endif
    }   // for (int i = 0; i < N; i++)

  for (int i = 0; i < (int)kdq_size(Q); i++)
    { elem = kdq_at(Q, i);
      cprofile[elem.seg_id].nw = plen-elem.b;   // TODO: 両端 WSIZE bp では他のリードで count MM がどう選ばれるか分からない。なので、Q に残っているものはすべて seed にする？
#ifdef DEBUG_QUEUE
      fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Last WSIZE bp)\n",
                      elem.seg_id,elem.cnt,elem.b,cprofile[elem.seg_id].nw);
#endif
    }

#ifdef DEBUG_SEGMENT
  fprintf(stderr,"# of segments = %d\n",N);
  for (int i = 0; i < N; i++)
    { seg = cprofile[i];
      fprintf(stderr,"seg[%d] (%d, %d) cnt = %d, # = %d\n",i,
                     seg.b,seg.e,seg.cnt,seg.nw);
    }
#endif

#ifdef TIME_SEED
  fprintf(stderr,"cm count\n");
  timeTo(stderr,false);
#endif

  // Pick up segments as "seed segments" (i.e. segments with the same count from each of which
  // minimizer(s) are selected) until the read gets covered by either invalid segments (i.e.
  // non-`C`-segments) or seed segments.
  M = 0;
  for (int i = 0; i < N; i++)
    { seg = cprofile[i];
      if (seg.cnt == -1)
        { mintvl[M].b = seg.b;
          mintvl[M].e = seg.e;
          M++;
        }
    }

#ifdef DEBUG_SEGMENT
  fprintf(stderr,"# of invalid segments = %d\n",M);
  for (int i = 0; i < M; i++)
    fprintf(stderr,"mintvl[%d] (%d, %d)\n",i,mintvl[i].b,mintvl[i].e);
#endif

  if (mintvl[0].b == 0 && mintvl[0].e == plen)
    goto seed_exit;

  // Process from segment with the largest # of windows until the masked interval contains [0..plen]
  qsort(cprofile,N,sizeof(seg_t),compare_cprofile);

#ifdef DEBUG_SEGMENT
  fprintf(stderr,"Sorted segments:\n");
  for (int i = 0; i < N; i++)
    fprintf(stderr,"seg[%d] = [(%d, %d) cnt = %d, # = %d]\n",
                   i,cprofile[i].b,cprofile[i].e,cprofile[i].cnt,cprofile[i].nw);
  fprintf(stderr,"\n");
#endif

  for (int i = 0; i < N; i++)
    { seg = cprofile[i];
#ifdef DEBUG_SELECT
      fprintf(stderr,"i = %d, seg = (%d, %d), cnt = %d, # = %d\n",
                     i,seg.b,seg.e,seg.cnt,seg.nw);
#endif
      if (seg.nw == WSIZE || !is_contained(mintvl,M,seg.b,seg.e))
        { M = add_intvl(mintvl,M,MAX(0,seg.b-WSIZE),MIN(seg.e+WSIZE,plen));
          cprofile[i].is_seed = true;

          // Choose only minimizer(s)
          int min_hash = MOD;
          for (int j = seg.b; j < seg.e; j++)
            min_hash = MIN(hash[j],min_hash);
          for (int j = seg.b; j < seg.e; j++)
            if (hash[j] == min_hash)
              sasgn[j] = -2;

#ifdef DEBUG_SELECT
          fprintf(stderr,"-> Seed [(%d, %d) cnt = %d, nw = %d]\n",
                         seg.b,seg.e,seg.cnt,seg.nw);
#endif

          if (mintvl[0].b == 0 && mintvl[0].e == plen)
            break;
        }
    }

#ifdef DEBUG_SELECT
  int n = 0;
  for (int i = 0; i < N; i++)
    if (cprofile[i].is_seed)
      n++;
  fprintf(stderr,"\n# of %c-seed segments = %d\n",C,n);
  for (int i = 0; i < N; i++)
    { seg = cprofile[i];
      if (seg.is_seed)
        fprintf(stderr,"Seed seg[%d] = [(%d, %d) cnt = %d, nw = %d]\n",
                       i,seg.b,seg.e,seg.cnt,seg.nw);
    }
  if (!(mintvl[0].b == 0 && mintvl[0].e == plen))
    { fprintf(stderr,"Invalid segments | (Seed segments +- WSIZE) do not cover the read\n");
      exit(1);
    }
#endif

#ifdef TIME_SEED
  fprintf(stderr,"cm select\n");
  timeTo(stderr,false);
#endif

seed_exit:
  // Empty the queue
  kdq_size(Q) = 0;

  return;
}

/* NOTE: len(profile) == len(class) == plen
         len(_seq) == plen + K - 1
   NOTE: class[i] in {'E', 'H', 'D', 'R'}
*/
void find_seeds(kdq_t(hmer_t) *Q,
                const char    *_seq,
                const char    *class,
                const uint16  *profile,
                seg_t         *cprofile,
                int           *hash,
                int           *sasgn,
                intvl_t       *mintvl,
                const int      plen,
                const int      K)
{
#if defined(DEBUG_SEED) || defined(DEBUG_QUEUE) || defined(INFO_SEED)
  fprintf(stderr,"\n");
#endif

  const char *seq = _seq+K-1;   // kmer @ i on profile == seq[i-K+1]..seq[i]

#ifdef TIME_SEED
  fprintf(stderr,"seed init\n");
  timeTo(stderr,false);
#endif

  // Compute canonical hash for every k-mer
  kmer_hash(seq,hash,plen,K);

#ifdef TIME_SEED
  fprintf(stderr,"hash\n");
  timeTo(stderr,false);
#endif

  // Find seeds from first H-mers and then D-mers
  _find_seeds(Q,profile,seq,class,hash,sasgn,cprofile,mintvl,plen,K,'H');
  _find_seeds(Q,profile,seq,class,hash,sasgn,cprofile,mintvl,plen,K,'D');   // TODO: density adjustment based on the H-MMs

  // change flag value into seed info
  for (int i = 0; i < plen; i++)
    { if (sasgn[i] == -2)         // H,D-seed
        sasgn[i] = class[i];
      else if (class[i] == 'E')   // Sequencing error
        sasgn[i] = 'E';
      else                        // Others (Non-seed H,D and R)
        sasgn[i] = 'R';
    }

#ifdef DEBUG_SEED
  for (int i = 0; i < plen; i++)
    if (sasgn[i] == 'H' || sasgn[i] == 'D')
      fprintf(stderr,"seed(%c) @ %5d: kmer = %.*s, count = %d\n",
                     class[i],i,K,seq+i-K+1,profile[i]);
#endif

#ifdef DEBUG_ITER
  int c = 0, e = 0;
  for (int i = 0; i < plen; i++)
    { if (sasgn[i] == 'H' || sasgn[i] == 'D')
        c++;
      else if (sasgn[i] == 'E')
        e++;
    }
  fprintf(stderr,", %3d seeds, %3d errors",c,e);
#endif

  return;
}
