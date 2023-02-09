#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "nthash.h"
#include "ClassPro.h"

// #undef DEBUG_HASH
// #define DEBUG_COMPRESS
// #define DEBUG_QUEUE
// #define DEBUG_SEED
// #define DEBUG_SEGMENT
// #define DEBUG_SELECT

// #define DEBUG_REP
// #define DEBUG_COMPRESS_REP
// #define DEBUG_QUEUE_REP
// #define DEBUG_SEGMENT_REP
// #define DEBUG_SELECT_REP

// #undef TIME_SEED

#define WSIZE 1000
#define WSIZE_REP 200
#define BOUNDARY_UNIQ_LEN 2000
#define MOD 2147483647

static void kmer_hash(const char *seq,
                      int        *hash,
                      int         plen,
                      int         K)
{ uint64_t hVal, fhVal, rhVal;   // canonical, forward, and reverse-strand hash values

  fhVal = rhVal = 0;
  hVal = NTC64_b(seq-K+1,K,&fhVal,&rhVal);   // initial hash value
  hash[0] = (int)(hVal % MOD);   // TODO: use int64 as-is?

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

/* Compress consecutive tie counts in `profile` into a single segment.
   Treat only k-mers with class `C`.
   Do not treat k-mers at position i where sasgn[i] == -10 (i.e. highly repetitive region)
*/
static int compress_profile(const uint16 *profile,
                            const char   *class,
                            const int    *sasgn,
                            seg_t        *cprofile,
                            const int     plen,
                            const char    C)
{ int N = 0;
  int b = 0, e = 1;
  // bool prev_valid = (class[0] == C && sasgn[0] > -10) ? true : false;
  bool prev_valid = (class[0] == C) ? true : false;   // FIXME: dev mod

  while (e < plen)
    { if (!prev_valid)
        { // while (e < plen && (class[e] != C || sasgn[e] <= -10))
          while (e < plen && class[e] != C)   // FIXME: dev mod
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
          // prev_valid = (class[b] == C && sasgn[b] > -10) ? true : false;
          prev_valid = (class[b] == C) ? true : false;   // FIXME: dev mod
        }
    }

#ifdef DEBUG_COMPRESS
  fprintf(stderr,"[compress_profile, %c] # of segments = %d\n",C,N);
  for (int i = 0; i < N; i++)
    fprintf(stderr,"[compress_profile, %c] %d: (%d, %d) cnt = %d\n",
                   C,i,cprofile[i].b,cprofile[i].e,cprofile[i].cnt);
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

static void _find_seeds(kdq_t(hmer_t) *Q,
                        const uint16  *profile,
                        const char    *class,
                        const int     *hash,
                        int           *sasgn,
                        seg_t         *cprofile,
                        intvl_t       *mintvl,
                        const int      plen,
                        const char     C)
{ int N, M;
  hmer_t now, first, elem;
  seg_t seg;

  N = compress_profile(profile,class,sasgn,cprofile,plen,C);

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
        { // Add a new element into the deque
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
                          fprintf(stderr,"seg[%d] (cnt=%d,b=%d) nw <- %d (Found larger count %d @ b=%d)\n",
                                         elem.seg_id,elem.cnt,elem.b,cprofile[elem.seg_id].nw,now.cnt,now.b);
#endif
                        }
                      else
                        { cprofile[elem.seg_id].nw = elem.cnt;   // TODO: this is ad-hoc. better 2nd round seeding?
                        }
                    }
                  kdq_size(Q) = 0;
                }
            }
          while (kdq_size(Q) > 0)
            { elem = kdq_at(Q, kdq_size(Q) - 1);
              if (elem.cnt < now.cnt)
                { cprofile[elem.seg_id].nw = elem.cnt;   // TODO: this is ad-hoc.
                  kdq_size(Q)--;
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

#ifdef DEBUG_QUEUE
      if (kdq_size(Q) == 0) continue;
      fprintf(stderr,"@ %5d (%c): |Q| = %ld, Q =",i,class[i],kdq_size(Q));
      for (int j = 0; j < (int)kdq_size(Q); j++) {
        fprintf(stderr,"  %d @ seg[%d] (b=%d)",
                       kdq_at(Q, j).cnt,kdq_at(Q, j).seg_id,kdq_at(Q, j).b);
      }
      fprintf(stderr,"\n");
#endif
    }   // for (int i = 0; i < N; i++)

    // Treat all the segments remaining in the queue as OoR
    while (kdq_size(Q) > 0)
      { first = kdq_at(Q, 0);
        if (last_oor)
          { cprofile[first.seg_id].nw = MIN(first.b-last_oor_pos+1,WSIZE);
#ifdef DEBUG_QUEUE
        fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Boundary OoR (last_oor) @ %d)\n",
                        first.seg_id,first.cnt,first.b,cprofile[first.seg_id].nw,last_oor_pos);
#endif
          }
        else
          { cprofile[first.seg_id].nw = WSIZE;
#ifdef DEBUG_QUEUE
        fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Boundary OoR)\n",
                        first.seg_id,first.cnt,first.b,cprofile[first.seg_id].nw);
#endif
          }
        if (kdq_size(Q) > 1 && first.cnt > kdq_at(Q, 1).cnt)
          { last_oor_pos = first.e;
          }
        kdq_shift(hmer_t, Q);
        last_oor = true;
      }

#ifdef DEBUG_SEGMENT
  // Show all segments with # of MM windows information
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

  // "Invalid" segments
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

  if (M > 0 && mintvl[0].b == 0 && mintvl[0].e == plen)
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

  { int i = 0;
    for (; i < N; i++)
      { seg = cprofile[i];
        if (seg.nw < WSIZE)
          break;
#ifdef DEBUG_SELECT
        fprintf(stderr,"i = %d, seg = (%d, %d), cnt = %d, # = %d\n",
                      i,seg.b,seg.e,seg.cnt,seg.nw);
#endif
        M = add_intvl(mintvl,M,MAX(0,seg.b-WSIZE),MIN(seg.e+WSIZE,plen));
        cprofile[i].is_seed = true;

        // Choose only minimizer(s)
        int min_hash = MOD;
        for (int j = seg.b; j < seg.e; j++)
          min_hash = MIN(hash[j],min_hash);
        for (int j = seg.b; j < seg.e; j++)
          if (hash[j] == min_hash)
            { sasgn[j] = -2;
#ifdef DEBUG_SEED
              fprintf(stderr,"seed(%c) @ %5d: count = %d, nw = %d\n",
                               C,j,profile[j],seg.nw);
#endif
            }

#ifdef DEBUG_SELECT
        fprintf(stderr,"-> Seed [(%d, %d) cnt = %d, nw = %d]\n",
                         seg.b,seg.e,seg.cnt,seg.nw);
#endif
      }
    while (i < N)
      {
#ifdef DEBUG_SELECT
        seg = cprofile[i];
        fprintf(stderr,"i = %d, seg = (%d, %d), cnt = %d, # = %d\n",
                      i,seg.b,seg.e,seg.cnt,seg.nw);
#endif
        int ii;
        for (ii = i; ii < N && cprofile[i].nw == cprofile[ii].nw; ii++)
          { if (!is_contained(mintvl,M,cprofile[ii].b,cprofile[ii].e))
              cprofile[ii].is_seed = true;
          }
        for (ii = i; ii < N && cprofile[i].nw == cprofile[ii].nw; ii++)
          { seg = cprofile[ii];
            if (seg.is_seed)
              { M = add_intvl(mintvl,M,MAX(0,seg.b-WSIZE),MIN(seg.e+WSIZE,plen));
                // cprofile[i].is_seed = true;
                // TODO: pick up other segments with the same count as well

                // Choose only minimizer(s)
                int min_hash = MOD;
                for (int j = seg.b; j < seg.e; j++)
                  min_hash = MIN(hash[j],min_hash);
                for (int j = seg.b; j < seg.e; j++)
                  if (hash[j] == min_hash)
                    { sasgn[j] = -2;
#ifdef DEBUG_SEED
                      fprintf(stderr,"seed(%c) @ %5d: count = %d, nw = %d\n",
                                    C,j,profile[j],seg.nw);
#endif
                    }

#ifdef DEBUG_SELECT
                fprintf(stderr,"-> Seed [(%d, %d) cnt = %d, nw = %d]\n",
                              seg.b,seg.e,seg.cnt,seg.nw);
#endif

              }
          }
        if (M > 0 && mintvl[0].b == 0 && mintvl[0].e == plen)
          break;
        i = ii;
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

  if (!(mintvl[0].b == 0 && mintvl[0].e == plen))   // Should never enter here
    { fprintf(stderr,"[ERROR] Invalid segments | (Seed segments +- WSIZE) do not cover the read\n");
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

/* Find repeat intervals on a profile.
   Unique intervals := consecutive H/D-mers longer than K*2.5 bp.
   Repeat intervals := [0, plen) - (unique intervals)
*/
static void anno_repeat(int *sasgn,
                        const char *class,
                        const int plen,
                        const int K,
                        FILE *ranno,
                        FILE *rdata,
                        int64 *ridx)
{ const int MIN_UNIQ_LEN = K * 2.5;

  for (int i = 0; i < plen; i++)
    sasgn[i] = -10;

  // Rescue non-R regions longer than `MIN_UNIQ_LEN`, while ignoring E-only intervals.
  // NOTE: Do not forget to add Km1 when output `b,e` to DAZZ_TRACK.
  int b = 0, e;
  bool in_R = (class[0] == 'R') ? true : false;
  int n_normal = (class[0] == 'H' || class[0] == 'D') ? 1 : 0;
  for (e = 1; e < plen; e++)
    { if (in_R)
        { if (class[e] != 'R')
            { b = e;
              in_R = false;
              n_normal = (class[e] == 'H' || class[e] == 'D') ? 1 : 0;
            }
        }
      else
        { if (class[e] == 'R')
            { if (n_normal >= MIN_UNIQ_LEN)
                { for (int i = b; i < e; i++)
                    sasgn[i] = 0;
#ifdef DEBUG_REP
                  fprintf(stderr,"[anno_repeat] Recovering [%d..%d) (>MIN_UNIQ_LEN)\n",b,e);
#endif
                }
              in_R = true;
            }
          else if (class[e] == 'H' || class[e] == 'D')
            n_normal++;
        }
    }
  if (!in_R)
    { if (n_normal >= MIN_UNIQ_LEN)
        { for (int i = b; i < e; i++)
            sasgn[i] = 0;
#ifdef DEBUG_REP
          fprintf(stderr,"[anno_repeat] Recovering [%d..%d) (>MIN_UNIQ_LEN)\n",b,e);
#endif
        }
    }
  // Up to here `b` and `e` are defined on [0..plen)

  // Write to DAZZ_TRACK
  int n_rintvl = 0;
  in_R = (sasgn[0] == -10) ? true : false;
  b = K-1;
  for (int i = 1; i < plen; i++)
    { if (!in_R)
        { if (sasgn[i] == -10)
            { b = i+K-1;
              in_R = true;
            }
        }
      if (in_R)
        { if (sasgn[i] != -10)
            { e = i+K-1;
#ifndef DEBUG_SINGLE
              fwrite(&b,sizeof(int),1,rdata);
              fwrite(&e,sizeof(int),1,rdata);
#endif
              n_rintvl++;
              in_R = false;
#ifdef DEBUG_REP
              fprintf(stderr,"[anno_repeat] R-intvl: [%d .. %d)\n",b,e);
#endif
            }
        }
    }
  if (in_R)
    { e = plen+K-1;
#ifndef DEBUG_SINGLE
      fwrite(&b,sizeof(int),1,rdata);
      fwrite(&e,sizeof(int),1,rdata);
#endif
      n_rintvl++;
#ifdef DEBUG_REP
      fprintf(stderr,"[anno_repeat] R-intvl: [%d .. %d)\n",b,e);
#endif
    }
#ifndef DEBUG_SINGLE
  *ridx += (n_rintvl*2*sizeof(int));
  fwrite(ridx,sizeof(int64),1,ranno);
#endif

  // Annotate non-boundary repeat/errors as -11   // TODO: no need of this anymore?
  int l = BOUNDARY_UNIQ_LEN;
  while (l < plen && sasgn[l] == -10) l++;
  int r = plen-BOUNDARY_UNIQ_LEN;
  while (r >= 0 && sasgn[r] == -10) r--;
#ifdef DEBUG_REP
  fprintf(stderr,"[anno_repeat] non-boundary interval = [%d, %d]\n",l,r);
#endif
  for (int i = l; i < r; i++)
    if (sasgn[i] == -10)
      sasgn[i] = -11;

#ifdef DEBUG_REP
  fprintf(stderr,"[anno_repeat] # of R-intvals = %d\n",n_rintvl);
#endif

  return;
}


/* Compress consecutive tie counts in `profile` into a single segment.
   Treat only non-Error k-mers.
   Treat only k-mers at position i where sasgn[i] == -10 (i.e. highly repetitive region)
*/
static int compress_profile_rep(const uint16 *profile,
                                const char   *class,
                                const int    *sasgn,
                                seg_t        *cprofile,
                                const int     plen)
{ int N = 0;
  int b = 0, e = 1;
  bool prev_valid = (sasgn[0] <= -10 && class[0] != 'E') ? true : false;

// #ifdef DEBUG_COMPRESS_REP
//   fprintf(stderr,"@ %5d: sasgn = %3d, class = %c\n",0,sasgn[0],class[0]);
// #endif

  while (e < plen)
    {
#ifdef DEBUG_COMPRESS_REP
      fprintf(stderr,"@ %5d: sasgn = %3d, class = %c\n",e,sasgn[e],class[e]);
#endif
      if (!prev_valid)
        { while (e < plen && (sasgn[e] > -10 || class[e] == 'E'))
            {
#ifdef DEBUG_COMPRESS_REP
              fprintf(stderr,"@ %5d: sasgn = %3d, class = %c\n",e,sasgn[e],class[e]);
#endif
              e++;
            }
#ifdef DEBUG_COMPRESS_REP
          fprintf(stderr,"[compress_profile_rep] Invalid segment: [%d..%d)\n",b,e);
#endif
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
#ifdef DEBUG_COMPRESS_REP
          fprintf(stderr,"[compress_profile_rep] Valid segment: [%d..%d) cnt = %5d\n",
                         b,e,profile[e-1]);
#endif
          cprofile[N].b = b;
          cprofile[N].e = e;
          cprofile[N].cnt = profile[e-1];
          cprofile[N].nw = 0;
          cprofile[N].is_seed = false;
          N++;
          b = e;
          e++;
          prev_valid = (sasgn[b] <= -10 && class[b] != 'E') ? true : false;
        }
    }

#ifdef DEBUG_COMPRESS_REP
  fprintf(stderr,"[compress_profile_rep] # of segments = %d\n",N);
  for (int i = 0; i < N; i++)
    fprintf(stderr,"[compress_profile_rep] %d: (%d, %d) cnt = %d\n",
                   i,cprofile[i].b,cprofile[i].e,cprofile[i].cnt);
#endif

  return N;
}

static void _find_seeds_rep(kdq_t(hmer_t) *Q,
                            const uint16  *profile,
                            const char    *class,
                            const int     *hash,
                            int           *sasgn,
                            seg_t         *cprofile,
                            intvl_t       *mintvl,
                            const int      plen)
{ int N, M;
  hmer_t now, first, elem;
  seg_t seg;

  N = compress_profile_rep(profile,class,sasgn,cprofile,plen);

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

  // Uniformly sparse count **minimizers** using sliding window algorithm
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
              if (first.cnt > now.cnt)   // all elements are wiped out from Q
                { last_oor = false;
                  for (int j = 0; j < (int)kdq_size(Q); j++)
                    { elem = kdq_at(Q, j);
                      if (first.cnt == elem.cnt)
                        { cprofile[elem.seg_id].nw = MIN(now.b - elem.b,WSIZE_REP);
#ifdef DEBUG_QUEUE_REP
                          fprintf(stderr,"seg[%d] (cnt=%d,b=%d) nw <- %d (Found smaller count)\n",
                                         elem.seg_id,elem.cnt,elem.b,cprofile[elem.seg_id].nw);
#endif
                        }
                      else
                        { // break;
                          cprofile[elem.seg_id].nw = MAX(WSIZE_REP-elem.cnt,0);   // TODO: this is very ad-hoc. Need a better way...
                        }
                    }
                  kdq_size(Q) = 0;
                }
            }
          while (kdq_size(Q) > 0)
            { elem = kdq_at(Q, kdq_size(Q) - 1);
              if (elem.cnt > now.cnt)
                { cprofile[elem.seg_id].nw = MAX(WSIZE_REP-elem.cnt,0);   // TODO: this is very ad-hoc. Need a better way...
                  kdq_size(Q)--;
                }
              else
                break;
            }
          kdq_push(hmer_t, Q, now);
        }

      if (kdq_size(Q) == 0) continue;

      // Remove out-of-range elements
      while (kdq_size(Q) > 0 && kdq_at(Q, 0).b <= seg.b - WSIZE_REP)
        { first = kdq_at(Q, 0);
          if (last_oor)
            { cprofile[first.seg_id].nw = MIN(first.b-last_oor_pos+1,WSIZE_REP);
#ifdef DEBUG_QUEUE_REP
          fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Out of range & last OoR @ %d)\n",
                          first.seg_id,first.cnt,first.b,cprofile[first.seg_id].nw,last_oor_pos);
#endif
            }
          else
            { cprofile[first.seg_id].nw = WSIZE_REP;
#ifdef DEBUG_QUEUE_REP
          fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Out of range)\n",
                          first.seg_id,first.cnt,first.b,cprofile[first.seg_id].nw);
#endif
            }
          if (kdq_size(Q) > 1 && first.cnt < kdq_at(Q, 1).cnt)
            { last_oor_pos = first.e;
            }
          kdq_shift(hmer_t, Q);
          last_oor = true;
        }

#ifdef DEBUG_QUEUE_REP
      if (kdq_size(Q) == 0) continue;
      fprintf(stderr,"@ %5d (%c): |Q| = %ld, Q =",i,class[i],kdq_size(Q));
      for (int j = 0; j < (int)kdq_size(Q); j++) {
        fprintf(stderr,"  %d @ seg[%d] (b=%d)",
                       kdq_at(Q, j).cnt,kdq_at(Q, j).seg_id,kdq_at(Q, j).b);
      }
      fprintf(stderr,"\n");
#endif
    }   // for (int i = 0; i < N; i++)

    // Treat all the segments remaining in the queue as OoR
    while (kdq_size(Q) > 0)
      { first = kdq_at(Q, 0);
        if (last_oor)
          { cprofile[first.seg_id].nw = MIN(first.b-last_oor_pos+1,WSIZE_REP);
#ifdef DEBUG_QUEUE
        fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Boundary OoR (last_oor) @ %d)\n",
                        first.seg_id,first.cnt,first.b,cprofile[first.seg_id].nw,last_oor_pos);
#endif
          }
        else
          { cprofile[first.seg_id].nw = WSIZE_REP;
#ifdef DEBUG_QUEUE
        fprintf(stderr,"seg[%d] (cnt=%d, b=%d) nw <- %d (Boundary OoR)\n",
                        first.seg_id,first.cnt,first.b,cprofile[first.seg_id].nw);
#endif
          }
        if (kdq_size(Q) > 1 && first.cnt > kdq_at(Q, 1).cnt)
          { last_oor_pos = first.e;
          }
        kdq_shift(hmer_t, Q);
        last_oor = true;
      }

#ifdef DEBUG_SEGMENT_REP
  // Show all segments with # of MM windows information
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

#ifdef DEBUG_SEGMENT_REP
  fprintf(stderr,"# of invalid segments = %d\n",M);
  for (int i = 0; i < M; i++)
    fprintf(stderr,"mintvl[%d] (%d, %d)\n",i,mintvl[i].b,mintvl[i].e);
#endif

  if (M > 0 && mintvl[0].b == 0 && mintvl[0].e == plen)
    goto seed_exit_rep;

  // Process from segment with the largest # of windows until the masked interval contains [0..plen]
  qsort(cprofile,N,sizeof(seg_t),compare_cprofile);

#ifdef DEBUG_SEGMENT_REP
  fprintf(stderr,"Sorted segments:\n");
  for (int i = 0; i < N; i++)
    fprintf(stderr,"seg[%d] = [(%d, %d) cnt = %d, # = %d]\n",
                   i,cprofile[i].b,cprofile[i].e,cprofile[i].cnt,cprofile[i].nw);
  fprintf(stderr,"\n");
#endif

  { int i = 0;
    for (; i < N; i++)
      { seg = cprofile[i];
        if (seg.nw < WSIZE_REP)
          break;
#ifdef DEBUG_SELECT_REP
        fprintf(stderr,"i = %d, seg = (%d, %d), cnt = %d, # = %d\n",
                       i,seg.b,seg.e,seg.cnt,seg.nw);
#endif
        M = add_intvl(mintvl,M,MAX(0,seg.b-WSIZE_REP),MIN(seg.e+WSIZE_REP,plen));
        cprofile[i].is_seed = true;

        // Choose only minimizer(s)
        int min_hash = MOD;
        for (int j = seg.b; j < seg.e; j++)
          min_hash = MIN(hash[j],min_hash);
        for (int j = seg.b; j < seg.e; j++)
          if (hash[j] == min_hash)
            { sasgn[j] = -3;
#ifdef DEBUG_SEED
              fprintf(stderr,"seed(R) @ %5d: count = %d, nw = %d\n",
                             j,profile[j],seg.nw);
#endif
            }

#ifdef DEBUG_SELECT_REP
        fprintf(stderr,"-> Seed [(%d, %d) cnt = %d, nw = %d]\n",
                        seg.b,seg.e,seg.cnt,seg.nw);
#endif
      }
    while (i < N)
      {
#ifdef DEBUG_SELECT_REP
        seg = cprofile[i];
        fprintf(stderr,"i = %d, seg = (%d, %d), cnt = %d, # = %d\n",
                       i,seg.b,seg.e,seg.cnt,seg.nw);
#endif
        int ii;
        for (ii = i; ii < N && cprofile[i].nw == cprofile[ii].nw; ii++)
          { if (!is_contained(mintvl,M,cprofile[ii].b,cprofile[ii].e))
              cprofile[ii].is_seed = true;
          }
        for (ii = i; ii < N && cprofile[i].nw == cprofile[ii].nw; ii++)
          { seg = cprofile[ii];
            if (seg.is_seed)
              { M = add_intvl(mintvl,M,MAX(0,seg.b-WSIZE_REP),MIN(seg.e+WSIZE_REP,plen));
                // cprofile[i].is_seed = true;

                // Choose only minimizer(s)
                int min_hash = MOD;
                for (int j = seg.b; j < seg.e; j++)
                  min_hash = MIN(hash[j],min_hash);
                for (int j = seg.b; j < seg.e; j++)
                  if (hash[j] == min_hash)
                    { sasgn[j] = -3;
#ifdef DEBUG_SEED
                      fprintf(stderr,"seed(R) @ %5d: count = %d, nw = %d\n",
                                    j,profile[j],seg.nw);
#endif
                    }

#ifdef DEBUG_SELECT_REP
                fprintf(stderr,"-> Seed [(%d, %d) cnt = %d, nw = %d]\n",
                              seg.b,seg.e,seg.cnt,seg.nw);
#endif
              }
          }
        if (M > 0 && mintvl[0].b == 0 && mintvl[0].e == plen)
          break;
        i = ii;
      }
  }

#ifdef DEBUG_SELECT_REP
  int n = 0;
  for (int i = 0; i < N; i++)
    if (cprofile[i].is_seed)
      n++;
  fprintf(stderr,"\n# of %c-seed segments = %d\n",'R',n);
  for (int i = 0; i < N; i++)
    { seg = cprofile[i];
      if (seg.is_seed)
        fprintf(stderr,"Seed seg[%d] = [(%d, %d) cnt = %d, nw = %d]\n",
                       i,seg.b,seg.e,seg.cnt,seg.nw);
    }

  // TODO: This happens. Fix.
  // if (!(mintvl[0].b == 0 && mintvl[0].e == plen))   // Should never enter here
  //   { fprintf(stderr,"[ERROR] Invalid segments | (Seed segments +- WSIZE) do not cover the read\n");
  //     exit(1);
  //   }
#endif

#ifdef TIME_SEED
  fprintf(stderr,"cm select\n");
  timeTo(stderr,false);
#endif

seed_exit_rep:
  // Empty the queue
  kdq_size(Q) = 0;

  return;
}

/* NOTE: len(profile) == len(class) == plen
         len(_seq) == plen + K - 1

   NOTE: class[i] in {'E', 'H', 'D', 'R'}

   NOTE: Meaning of the value of `sasgn` (for each position (= k-mer)):   // TODO: change to bit operation with macros?
     -11 = non-boundary repeats/errors
     -10 = boundary repeats1/errors
      -2 = seed fixed elsewhere
      -1 = fixed seed
       0 = candidate
       1 = seed
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
                const int      K,
                FILE          *ranno,
                FILE          *rdata,
                int64         *ridx)
{
#if defined(DEBUG_ITER) || defined(DEBUG_SEED) || defined(DEBUG_QUEUE)
  fprintf(stderr,"\n");
#endif

  const char *seq = _seq+K-1;   // kmer @ i on profile == seq[i-K+1]..seq[i]

#ifdef TIME_SEED
  fprintf(stderr,"seed init\n");
  timeTo(stderr,false);
#endif

  // Annotate highly repetitive regions
  anno_repeat(sasgn,class,plen,K,ranno,rdata,ridx);

  // Compute canonical hash for every k-mer
  kmer_hash(seq,hash,plen,K);

#ifdef TIME_SEED
  fprintf(stderr,"hash\n");
  timeTo(stderr,false);
#endif

  // Find seeds using count maximizers and then sequence minimizers,
  // from H-mers and D-mers NOT in highly repetitive regions
  _find_seeds(Q,profile,class,hash,sasgn,cprofile,mintvl,plen,'H');
  _find_seeds(Q,profile,class,hash,sasgn,cprofile,mintvl,plen,'D');

  // Find seeds among non-Error, low-copy k-mers from highly repetitive regions
  _find_seeds_rep(Q,profile,class,hash,sasgn,cprofile,mintvl,plen);

  // change flag value into seed info
  for (int i = 0; i < plen; i++)
    { if (sasgn[i] == -2)        // H,D-seed in normal regions
        sasgn[i] = class[i];
      else if (sasgn[i] == -3)   // Seeds from highly repetitive regions
        sasgn[i] = 'R';
      else
        sasgn[i] = 'E';   // Do not use
    }

#ifdef DEBUG_ITER
  int n = 0, r = 0;
  for (int i = 0; i < plen; i++)
    { if (sasgn[i] == 'H' || sasgn[i] == 'D')
        n++;
      else if (sasgn[i] == 'R')
        r++;
    }
  fprintf(stderr,"%3d normal seeds, %3d repeat seeds\n",n,r);
#endif

  return;
}
