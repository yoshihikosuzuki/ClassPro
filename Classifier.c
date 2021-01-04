#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "libfastk.h"
#include "DB.h"

#define WRITE_ASCII_CLASS
#undef DEBUG_ITER
#undef DEBUG_PROB

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define P_ERROR 0.001
#define MAX_NITER 10
#define CTX_LEN 1000

static char *Usage = "<source_root> <diplo_depth<int>>";

char stoc[4] = {'E', 'R', 'H', 'D'};

double logfact[32768];

inline double logp_poisson(int x, int lambda) {
  return x * log((double)lambda) - lambda - logfact[x];
}

inline double logp_binom(int k, int n, double p) {
  return logfact[n] - logfact[k] - logfact[n-k] + k * log(p) + (n-k) * log(1-p);
}

void nn_intvl(int i, char s, int s_not, int N, char* asgn, int ret[2]) {
  int p = i - 1;
  if (s_not)
    while (p >= 0 && asgn[p] == s) p--;
  else
    while (p >= 0 && asgn[p] != s) p--;
  ret[0] = p;
  int n = i + 1;
  if (s_not)
    while (n < N && asgn[n] == s) n++;
  else
    while (n < N && asgn[n] != s) n++;
  ret[1] = n;
}

void estimate_true_counts(int i, char s, int N, int* chpt, char* asgn, uint16* profile, int ret[2]) {
  int l, csum;

  l = CTX_LEN;
  csum = 0;
  int p = i - 1;
  while (p >= 0 && l > 0) {
    if (asgn[p] == s)
      for (int j = chpt[p+1]-1; j >= chpt[p] && l > 0; j--, l--)
        csum += profile[j];
    p--;
  }
  ret[0] = l == CTX_LEN ? -1 : (int)((double)csum / (CTX_LEN - l));
  l = CTX_LEN;
  csum = 0;
  int n = i + 1;
  while (n < N && l > 0) {
    if (asgn[n] == s)
      for (int j = chpt[n]; j < chpt[n+1] && l > 0; j++, l--)
        csum += profile[j];
    n++;
  }
  ret[1] = l == CTX_LEN ? -1 : (int)((double)csum / (CTX_LEN - l));
}

double calc_logp_h(int i, uint16* profile, int* chpt, char* asgn, int N, int* depths) {
  int ib, ie, cib, cie;
  int nn_pos[2];
  int p, n;
  int est_cnt[2];
  int lc, rc;
  double logp_p, logp_l, logp_r;

  ib = chpt[i];
  ie = chpt[i + 1] - 1;
  cib = profile[ib];
  cie = profile[ie];
  // H < D
  nn_intvl(i, 3, 0, N, asgn, nn_pos);
  p = nn_pos[0]; n = nn_pos[1];
  if (p >= 0 && profile[chpt[p + 1] - 1] < cib)
    return -DBL_MAX;
  if (n < N && profile[chpt[n]] < cie)
    return -DBL_MAX;
  // TODO: interval length requirement if surrounded by D?
  // prior
  logp_p = logp_poisson(cib, depths[2]);
  // transition
  estimate_true_counts(i, 2, N, chpt, asgn, profile, est_cnt);
  lc = est_cnt[0]; rc = est_cnt[1];
  if (lc == -1 && rc == -1) {
    logp_l = logp_poisson(cib, depths[2]);
    logp_r = logp_poisson(cie, depths[2]);
    return logp_p + logp_l + logp_r;
  }
  // left transition
  if (lc >= 0)
    logp_l = logp_binom(MIN(lc, cib), MAX(lc, cib), 1 - P_ERROR);
  // right transition
  if (rc >= 0)
    logp_r = logp_binom(MIN(rc, cie), MAX(rc, cie), 1 - P_ERROR);
  if (lc == -1) logp_l = logp_r;
  if (rc == -1) logp_r = logp_l;
  return logp_l + logp_r + logp_p;
}

double calc_logp_d(int i, uint16* profile, int* chpt, char* asgn, int N, int* depths) {
  int ib, ie, cib, cie;
  int nn_pos[2];
  int p, n;
  int est_cnt[2];
  int lc, rc;
  double logp_p, logp_l, logp_r;

  ib = chpt[i];
  ie = chpt[i + 1] - 1;
  cib = profile[ib];
  cie = profile[ie];

  // H < D
  nn_intvl(i, 2, 0, N, asgn, nn_pos);
  p = nn_pos[0]; n = nn_pos[1];
  if (p >= 0 && profile[chpt[p + 1] - 1] > cib)
    return -DBL_MAX;
  if (n < N && profile[chpt[n]] > cie)
    return -DBL_MAX;
  // prior
  logp_p = logp_poisson(cib, depths[3]);
  // transition
  estimate_true_counts(i, 3, N, chpt, asgn, profile, est_cnt);
  lc = est_cnt[0]; rc = est_cnt[1];
  if (lc == -1 && rc == -1) {
    logp_l = logp_poisson(cib, depths[3]);
    logp_r = logp_poisson(cie, depths[3]);
    return logp_p + logp_l + logp_r;
  }
  // left transition
  if (lc >= 0)
    logp_l = logp_binom(MIN(lc, cib), MAX(lc, cib), 1 - P_ERROR);
  // right transition
  if (rc >= 0)
    logp_r = logp_binom(MIN(rc, cie), MAX(rc, cie), 1 - P_ERROR);
  if (lc == -1) logp_l = logp_r;
  if (rc == -1) logp_r = logp_l;
  return logp_l + logp_r + logp_p;
}

double calc_logp_r(int i, uint16* profile, int* chpt, char* asgn, int N, int* depths) {
  int ib, ie, cib, cie;
  int est_cnt[2];
  int lc, rc;
  double lrcov, rrcov;
  double logp_l, logp_r;

  ib = chpt[i];
  ie = chpt[i + 1] - 1;
  cib = profile[ib];
  cie = profile[ie];

  estimate_true_counts(i, 3, N, chpt, asgn, profile, est_cnt);
  lc = est_cnt[0]; rc = est_cnt[1];
  if (lc == -1 && rc == -1)
    lrcov = rrcov = (double)depths[3] * 1.25;
  else {
    if (lc == -1) lc = rc;
    if (rc == -1) rc = lc;
    lrcov = lc * 1.25; rrcov = rc * 1.25;
  }
  // emission
  if (cib >= lrcov || cie >= rrcov)
      return 0.;
  logp_l = logp_binom(cib, (int)lrcov, 1 - P_ERROR);
  logp_r = logp_binom(cie, (int)rrcov, 1 - P_ERROR);
  return logp_l + logp_r;
}

double calc_logp_e(int i, uint16* profile, int* chpt, char* asgn, int N, int* depths) {
  int ib, ie, cib, cie;
  int nn_pos[2];
  int p, n;
  int pe, nb, cpe, cnb;
  double logp_pl, logp_pr, logp_l, logp_r;

  ib = chpt[i];
  ie = chpt[i + 1] - 1;
  cib = profile[ib];
  cie = profile[ie];
  
  // TODO: interval length requirement?
  // prior
  logp_pl = logp_poisson(cib, depths[0]);
  logp_pr = logp_poisson(cie, depths[0]);
  // transition
  nn_intvl(i, 0, 1, N, asgn, nn_pos);
  p = nn_pos[0]; n = nn_pos[1];
  // left transition
  if (p < 0) {
    logp_l = 0;
  } else {
    pe = chpt[p + 1] - 1;
    cpe = profile[pe];
    if (cpe <= cib)
      return -DBL_MAX;
    logp_l = logp_binom(cib, cpe, P_ERROR);
  }
  // right transition
  if (n >= N) {
    logp_r = 0;
  } else {
    nb = chpt[n];
    cnb = profile[nb];
    if (cnb <= cie)
      return -DBL_MAX;
    logp_r = logp_binom(cie, cnb, P_ERROR);
  }
  return logp_l + logp_r + logp_pl + logp_pr;
}

int update_state(int i, uint16* profile, int* chpt, char* asgn, int N, int* depths) {
  char s, smax = 5;
  double logp, logpmax = -DBL_MAX;

  for (s = 0; s < 4; s++) {
    if (s == 0)
      logp = calc_logp_e(i, profile, chpt, asgn, N, depths);
    else if (s == 1)
      logp = calc_logp_r(i, profile, chpt, asgn, N, depths);
    else if (s == 2)
      logp = calc_logp_h(i, profile, chpt, asgn, N, depths);
    else
      logp = calc_logp_d(i, profile, chpt, asgn, N, depths);
#ifdef DEBUG_PROB
    printf("i=%d, s=%d, logp=%lf\n",i,s,logp);
#endif
    if (logp > logpmax) {
      smax = s;
      logpmax = logp;
    }
  }
  if (smax > 4) {
    printf("No valid probability for interval %d\n",i);
    exit(1);
  }
  if (smax != asgn[i]) {
#ifdef DEBUG_ITER
    printf("state updated @ %d: %d -> %d\n",i,asgn[i],smax);
#endif
    asgn[i] = smax;
    return 1;
  } else {
    return 0;
  }
}

int main(int argc, char *argv[])
{ Profile_Index *P;
  DAZZ_DB        _db, *db = &_db;
  FILE          *afile, *dfile;   // .class.anno, .class.data
  //FILE          *rafile, *rdfile;   // .krep.anno, .krep.data
  int            nreads;
  int            dcov;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("Classifier");

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    dcov = (int)strtol(argv[2],&eptr,10);
    if (argv[2] == eptr || *eptr != 0)
      { fprintf(stderr,"%s: diplo_depth must be an integer\n",Prog_Name);
        exit(1);
      }
  }

  //  Open DB or DAM, and if a DAM open also .hdr file
  { int status;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (db->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block\n",Prog_Name);
        exit (1);
      }
    Trim_DB(db);
    nreads = db->nreads;
  }

  //  Create class and krep tracks
  { int size;
    char   *path, *root;

    path   = PathTo(argv[1]);
    root  = Root(argv[1],".db");

    afile = Fopen(Catenate(path,"/.",root,".class.anno"),"w+");
    dfile = Fopen(Catenate(path,"/.",root,".class.data"),"w+");
    if (afile == NULL || dfile == NULL)
      exit (1);

    size = 8;
    fwrite(&nreads,sizeof(int),1,afile);
    fwrite(&size,sizeof(int),1,afile);

    /*rafile = Fopen(Catenate(path,"/.",root,".krep.anno"),"w+");
    rdfile = Fopen(Catenate(path,"/.",root,".krep.data"),"w+");
    if (rafile == NULL || rdfile == NULL)
      exit (1);

    size = 0;
    fwrite(&P->nreads,sizeof(int),1,rafile);
    fwrite(&size,sizeof(int),1,rafile);*/

    free(path);
    free(root);
  }

  // Load Profile for each read
  P = Open_Profiles(argv[1]);
  if (P == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  if (nreads != P->nreads)
    { fprintf(stderr,"# of reads must be same\n");
      exit(1);
    }

  { int     id;
    uint16 *profile;
    int     plen, tlen;
    int     km1;
    int64   idx/*, rdx*/;
    char   *track, *crack;   // classification tracks
    int    *chpt;
    char   *asgn;
    int     N;

    const int sigch = 3;
    int *depths = Malloc(4*sizeof(int),"Depths for initial assignments");
    // 0: E, 1: R, 2: H, 3: D
    depths[0] = 1;
    depths[1] = dcov * 2;
    depths[2] = dcov / 2;
    depths[3] = dcov;

    for (int n = 1; n < 32768; n++)
      logfact[n] = logfact[n-1] + log(n);

    km1 = P->kmer - 1;
    track = New_Read_Buffer(db);
    crack = track + km1;

    /*mask = (int *) Malloc(sizeof(int)*block->maxlen,"Allocating mask vector");
    cvec = (uint16 *) Malloc(sizeof(uint16)*block->maxlen,"Allocating mask vector");
    if (mask == NULL)
      exit (1);*/

    idx = 0;
    fwrite(&idx,sizeof(int64),1,afile);
    /*rdx = 0;
    fwrite(&rdx,sizeof(int64),1,rafile);*/

    plen    = 20000;
    profile = Malloc(plen*sizeof(uint16),"Profile array");
    chpt    = Malloc(plen*sizeof(int),"Change point array");
    asgn    = Malloc(plen*sizeof(char),"Assignment array");

    for (id = 1; id <= nreads; id++)
      { tlen = Fetch_Profile(P,(int64) id-1,plen,profile);
        if (tlen > plen)
          { plen    = 1.2*tlen + 1000;
            profile = Realloc(profile,plen*sizeof(uint16),"Profile array");
            chpt    = Realloc(chpt,plen*sizeof(int),"Change point array");
            asgn    = Realloc(asgn,plen*sizeof(char),"Assignment array");
            Fetch_Profile(P,(int64) id-1,plen,profile);
          }

#ifdef DEBUG_ITER
        printf("Read %d:\n",id);
#endif
        
        for (int i = 0; i < km1; i++)
          track[i] = 0;

        // split into smooth intervals
        N = 0;
        chpt[N++] = 0;
        for (int i = 1; i < tlen; i++)
          if (abs(MIN(profile[i],depths[1]) - MIN(profile[i-1],depths[1])) >= sigch)
            chpt[N++] = i;
        chpt[N] = tlen;
        // initial assignment
        for (int i = 0; i < N; i++) {
            int dmin = INT_MAX;
            // TODO: binary search is faster?
            for (int s = 0; s < 4; s++) {
                int d = abs(depths[s] - MAX(profile[chpt[i]], profile[chpt[i+1]-1]));
                if (d < dmin) {
                    asgn[i] = s;
                    dmin = d;
                }
            }
        }
#ifdef DEBUG_ITER
        for (int i = 0; i < N; i++)
          printf("intvl %d-%d: %d\n",chpt[i],chpt[i+1]-1,asgn[i]);
#endif

        // iteratively refine assignments
        int counter;
        for (counter = 0; counter < MAX_NITER; counter++) {
#ifdef DEBUG_ITER
          printf(".");
#endif
          int changed = 0;
          for (int i = 0; i < N; i++)
            if (update_state(i, profile, chpt, asgn, N, depths))
              changed = 1;
          if (!changed) break;
#ifdef DEBUG_ITER
          printf(".");
#endif
          changed = 0;
          for (int i = N - 1; i >= 0; i--)
            if (update_state(i, profile, chpt, asgn, N, depths))
              changed = 1;
          if (!changed) break;
        }
#ifdef DEBUG_ITER
        if (counter == MAX_NITER) printf("Not converged");
        printf("\n");
#endif

        // map to track
        for (int i = 0; i < N; i++) {
          for (int j = chpt[i]; j < chpt[i + 1]; j++) {
            crack[j] = asgn[i];
          }
        }

        // change assignments of high-copy R into E
        /*for (int i = 0; i < tlen; i++)
          if (crack[i] == 1 && profile[i] > depths[1])
            crack[i] = 0;*/

#ifdef WRITE_ASCII_CLASS
        printf("%d\t",id);
        for (int i = 0; i < tlen; i++)
          printf("%c",stoc[(unsigned char)crack[i]]);
        printf("\n");
#endif

        int l = tlen + km1;
        Compress_Read(l,track);
        int t = COMPRESSED_LEN(l);
        fwrite(track,1,t,dfile);
        idx += t;
        fwrite(&idx,sizeof(int64),1,afile);

        /*if (r >= 0)
          mask[n++] = r;
        fwrite(mask,sizeof(int),n,rdfile);
        rdx += n*sizeof(int);
        fwrite(&rdx,sizeof(int64),1,rafile);*/
      }
    //free(mask);
    free(track-1);
    free(profile);
  }

  Free_Profiles(P);

  //fclose(rdfile);
  //fclose(rafile);
  fclose(dfile);
  fclose(afile);

  Close_DB(db);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}