#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>

#include "libfastk.h"
#include "DB.h"
#include "ClassPro.h"

#define THREAD pthread_t

int       VERBOSE;
int       NTHREADS;
int       NREADS;
int       NPARTS;
int       FIND_SEEDS;
int       READ_LEN;

static void *kmer_class_thread(void *arg)
{ Class_Arg     *data   = (Class_Arg *)arg;
  Profile_Index *P      = data->P;
  Error_Model   *emodel = data->emodel;
  const int      K      = P->kmer;
  const int      Km1    = K-1;
  
  char        *track, *crack;         // Data for dfile (crack = track + km1)
  int64        idx;                   // Data for afile

  DAZZ_READ   *r;
  char        *seq;                   // Fetched sequence of the read
  uint16      *profile, *nprofile;    // Fetched profile of a read
  int          rlen, plen;            // `rlen` = `plen` + `Km1`

  Seq_Ctx     *ctx[N_WTYPE];          // Context lengths per position
  Seq_Ctx     *lctx, *_lctx, *rctx;

  P_Error     *perror, *cerror;       // Error probability per position

  double      *eta;                   // Parameter for PMM
  double       lambda[2];             // H-cov, D-cov

  Error_Intvl *eintvl[N_ETYPE];
  Intvl       *intvl;
  Rel_Intvl   *rintvl;
  int         *wall;                  // Wall positions
  char        *asgn;                  // Interval classifications
  int          N, M;                  // Number of walls/reliable intervals
  char        *buf;

  // .db
  DAZZ_DB    *db    = data->db;
  DAZZ_STUB  *stub  = data->stub;
  char      **flist = stub->prolog;
  int        *findx = stub->nreads;
  int         map   = 0;
  findx[-1] = 0;

  // Prepare buffers etc.
  { idx = 0;
    
    buf = Malloc((db->maxlen+1)*sizeof(char),"buf");
    for (int i = 0; i < Km1; i++)
      buf[i] = 'N';

    seq   = New_Read_Buffer(db);
    track = New_Read_Buffer(db);
    crack = track + Km1;
    for (int i = 0; i < Km1; i++)
      track[i] = 0;

    intvl    = Malloc(db->maxlen*sizeof(Intvl),"Interval array");
    rintvl   = Malloc(db->maxlen*sizeof(Rel_Intvl),"Reliable interval array");
    for (int i = SELF; i <= OTHERS; i++)
      eintvl[i] = Malloc(db->maxlen*sizeof(Error_Intvl),"Error intvl array");
    wall     = Malloc(db->maxlen*sizeof(int),"Wall array");
    asgn     = Malloc(db->maxlen*sizeof(char),"Interval assignment array");
    perror   = Malloc(db->maxlen*sizeof(P_Error),"Error prob");
    cerror   = Malloc(db->maxlen*sizeof(P_Error),"Error prob");
    eta      = Malloc(db->maxlen*2*sizeof(double),"PMM eta");
    profile  = Malloc(db->maxlen*sizeof(uint16),"Profile array");
    nprofile = Malloc(db->maxlen*sizeof(uint16),"Normal profile array");
    _lctx    = Malloc(db->maxlen*sizeof(Seq_Ctx),"Allocating left ctx vector");
    rctx     = Malloc(db->maxlen*sizeof(Seq_Ctx),"Allocating right ctx vector");

    lctx     = _lctx + Km1 - 1;
    _lctx[0][HP] = 1;
    _lctx[0][DS] = _lctx[0][TS] = _lctx[1][TS] = 0;
    ctx[DROP] = lctx;
    ctx[GAIN] = rctx;
  }

#ifdef DEBUG_SMALL
  for (int id = data->beg; id < data->beg+NREAD_SMALL; id++)
#else
  for (int id = data->beg; id < data->end; id++)
#endif
    { r = db->reads+id;
      rlen = r->rlen;

#if defined(DEBUG) || defined(DEBUG_CTX) || defined(DEBUG_ERROR) || defined(DEBUG_ITER)
      fprintf(stderr,"\nRead %d (%d bp): ",id+1,rlen);
      fflush(stderr);
#endif

      Load_Read(db,id,seq,2);

#ifdef DEBUG
      const int slen = strlen(seq);
      if (rlen != slen)
        { fprintf(stderr,"rlen (%d) != strlen(seq) (%d)\n",rlen,slen);
          exit(1);
        }
      if (rlen > db->maxlen)
        { fprintf(stderr,"rlen (%d) > db->maxlen (%d)\n",rlen,db->maxlen);
          exit(1);
        }
#endif

      calc_seq_context(_lctx,rctx,seq,rlen);

      plen = Fetch_Profile(P,(int64)id,db->maxlen,profile);

      if (rlen != plen+Km1)
        { fprintf(stderr,"Read %d: rlen (%d) != plen+Km1 (%d)\n",id+1,rlen,plen+Km1);
          exit(1);
        }

      find_wall(profile,plen,ctx,emodel,perror,cerror,eintvl,intvl,rintvl,wall,asgn,K,&N,&M);

#ifdef DEBUG_ITER
      int rtot = 0;
      for (int i = 0; i < M; i++)
        rtot += rintvl[i].j - rintvl[i].i;
      fprintf(stderr,"%d (/%d; %d %%) rel intvls, ",M,N,(int)(100.*rtot/plen));
#endif

      int nnorm = pmm_vi(profile,nprofile,plen,eta,lambda);

#ifdef DEBUG_ITER
      fprintf(stderr,"(H,D)=(%.lf,%.lf) (%d %% normal)\n",lambda[0],lambda[1],(int)(100.*nnorm/plen));
#endif

      classify_reliable(rintvl,M,intvl,N,plen,perror,cerror,(int)lambda[0],(int)lambda[1]);
      classify_unreliable(profile,plen,intvl,N,perror,(int)lambda[0],(int)lambda[1]);

      for (int i = 0; i < N; i++)
        for (int j = intvl[i].i; j < intvl[i].j; j++)
          crack[j] = intvl[i].asgn;

      remove_slip(profile,plen,ctx,crack);

/*#ifdef DEBUG_ITER
      fprintf(stderr,"  Final: ");
      for (int i = 0; i < plen; i++)
        fprintf(stderr,"%c",stoc[(unsigned char)crack[i]]);
      fflush(stderr);
#endif*/

      // .db
      while (id < findx[map-1])
        map -= 1;
      while (id >= findx[map])
        map += 1;

      int bufidx = Km1;
      for (int i = 0; i < plen; i++)
        buf[bufidx++] = stoc[(unsigned char)crack[i]];
      buf[bufidx] = '\0';

      fprintf(data->cfile,"@%s/%d/%d_%d\n%s\n+\n%s\n",flist[map],r->origin,r->fpulse,r->fpulse+r->rlen,seq,buf);

      // Output binary to DAZZ track
      int l = plen + Km1;
      Compress_Read(l,track);
      int t = COMPRESSED_LEN(l);
      fwrite(track,1,t,data->dfile);
      idx += t;
      fwrite(&idx,sizeof(int64),1,data->afile);
    }

#ifdef DEBUG_MERGE
      fprintf(stderr,"last idx=%lld @thread %d\n",idx,data->wch);
      fflush(stderr);
#endif

  // .db
  //if (strcmp(ext,".db") == 0)
    { Close_DB(db);
      Free_DB_Stub(stub);
    }

  fclose(data->cfile);
  fclose(data->afile);
  fclose(data->dfile);

  free(track-1);
  free(seq-1);
  free(profile);
  free(nprofile);
  free(_lctx);
  free(rctx);
  free(perror);
  free(cerror);
  for (int i = SELF; i <= OTHERS; i++)
    free(eintvl[i]);
  free(intvl);
  free(rintvl);
  free(wall);
  free(asgn);
  free(eta);
  free(buf);

  return (NULL);
}

int main(int argc, char *argv[])
{ char          *path, *root, *ext;
  char          *FK_ROOT;
  char         **fnames;
  int            nfiles;

  Profile_Index *P;
  Error_Model   *emodel;

  THREAD        *threads;
  Class_Arg     *paramc;
  Merge_Arg     *paramm;

  int            COVERAGE;

  // Parse options
  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("ClassPro");

    NTHREADS = 4;
    COVERAGE = -1;
    READ_LEN = DEFAULT_RLEN;
    FK_ROOT   = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vs")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
          case 'c':
            ARG_NON_NEGATIVE(COVERAGE,"Estimated k-mer coverage")
            break;
          case 'r':
            ARG_POSITIVE(READ_LEN,"Average read length")
            break;
          case 'N':
            FK_ROOT = argv[i]+2;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE    = flags['v'];
    FIND_SEEDS = flags['s'];

    if (argc < 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
    nfiles = argc-1;
    fnames = argv+1;

    int fid, idx;

    path = PathTo(fnames[0]);
    for (idx = 0; idx < N_EXT; idx++)
      { root  = Root(fnames[0],EXT[idx]);
        fid   = open(Catenate(path,"/",root,EXT[idx]),O_RDONLY);
        if (fid >= 0) break;
        free(root);
      }
    if (idx == N_EXT || fid < 0)
      { fprintf(stderr,"Cannot open %s as a .f{ast}[aq][.gz]|db file\n",fnames[0]);
        exit(1);
      }
    close(fid);
    ext = EXT[idx];

    if (FK_ROOT == NULL)
      FK_ROOT = Strdup(Catenate(path,"/",root,""),"FK_ROOT");

    if (VERBOSE)
      { fprintf(stderr,"# of input files = %d (first file: path = %s  root = %s  ext = %s)\n",nfiles,path,root,ext);
        fprintf(stderr,"FASTK root = %s\n",FK_ROOT);
      }
  }

  // Precompute
  { if ((P = Open_Profiles(FK_ROOT)) == NULL)
      { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,FK_ROOT);
        exit (1);
      }
    NREADS = P->nreads;
    NPARTS = (NREADS / NTHREADS) + (NREADS % NTHREADS == 0 ? 0 : 1);
    if (VERBOSE)
      fprintf(stderr,"NREADS=%d, NPARTS=%d\n",NREADS,NPARTS);

    precompute_probs();
    process_global_hist(FK_ROOT,COVERAGE);
    emodel = calc_init_thres();
  }

  // Invoke classification 
  { threads = Malloc(sizeof(THREAD)*NTHREADS,"Allocating class threads");
    paramc  = Malloc(sizeof(Class_Arg)*NTHREADS,"Allocating class args");
    paramm  = Malloc(sizeof(Merge_Arg)*N_OTYPE,"Allocate merge args");

    if (strcmp(ext,".db") == 0)
      { if (nfiles != 1)
          { fprintf(stderr,"# of input .db files must be 1\n");
            exit(1);
          }
        prepare_db(path,root,paramc);
      }
    else
      prepare_fx(fnames,nfiles,paramc,path,root);

    // NOTE: 14 for "/." and ".class.[anno|data]." and 10 for "`t`" and '\0'
    const int fnlen = strlen(path)+strlen(root)+14+10;
    for (int c = 0; c < N_OTYPE; c++)
      { paramm[c].fnames = Malloc(sizeof(char *)*NTHREADS,"Allocating fnames");
        for (int t = 0; t < NTHREADS; t++)
          { paramm[c].fnames[t] = Malloc(sizeof(char)*fnlen,"Allocating fname");
            sprintf(paramm[c].fnames[t],"%s%s%s%s.%d",path,osep[c],root,osuf[c],t+1);
          }
        paramm[c].final = Malloc(sizeof(char)*fnlen,"Allocating fname");
        sprintf(paramm[c].final,"%s%s%s%s",path,osep[c],root,osuf[c]);
        paramm[c].N         = NTHREADS;
        paramm[c].is_binary = obin[c];
      }

    for (int t = 0; t < NTHREADS; t++)
      { 
#ifndef DUP_PROFILE
        paramc[t].P = (t == 0) ? P : Clone_Profiles(P);
#else
        paramc[t].P = (t == 0) ? P : Open_Profiles(FK_ROOT);
#endif
        if (paramc[t].P == NULL)
          { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,FK_ROOT);
            exit (1);
          }

        paramc[t].emodel = emodel;
        paramc[t].beg    = t*NPARTS;
        paramc[t].end    = MIN((t+1)*NPARTS,NREADS);

        paramc[t].afile = Fopen(paramm[ANNO].fnames[t],"wb");
        paramc[t].dfile = Fopen(paramm[DATA].fnames[t],"wb");
        if (paramc[t].afile == NULL || paramc[t].dfile == NULL)
          { fprintf(stderr,"Cannot open .class.*.%d\n",t+1);
            exit (1);
          }

#ifndef PARALLEL_WRITE
        paramc[t].cfile  = Fopen(paramm[CLASS].fnames[t],"w");
        if (paramc[t].cfile == NULL)
          { fprintf(stderr,"Cannot open *.class.%d\n",t+1);
            exit (1);
          }
#endif
      }

    { const int idx  = 0;
      const int size = 8;
      fwrite(&NREADS,sizeof(int),1,paramc[0].afile);
      fwrite(&size,sizeof(int),1,paramc[0].afile);
      fwrite(&idx,sizeof(int64),1,paramc[0].afile);
    }

    if (VERBOSE)
      fprintf(stderr,"Classifying %d-mers%s...\n",P->kmer,FIND_SEEDS ? "" : " & Finding seeds");

    for (int t = 1; t < NTHREADS; t++)
      pthread_create(threads+t,NULL,kmer_class_thread,paramc+t);
    kmer_class_thread(paramc);
    for (int t = 1; t < NTHREADS; t++)
      pthread_join(threads[t],NULL);

    { int i;

#ifndef PARALLEL_WRITE
      const int c = CLASS;
#else
      const int c = DATA;
#endif

      if (VERBOSE)
        fprintf(stderr,"\nMerging files...\n");

      for (i = 0; i+1 < NTHREADS && c+i < ANNO; i++)
        pthread_create(threads+i+1,NULL,merge_files,paramm+c+i);
      for (; c+i < ANNO; i++)
        merge_files(paramm+c+i);
      merge_anno(paramm+ANNO);
      for (i = 0; i+1 < NTHREADS && c+i < ANNO; i++)
        pthread_join(threads[i+1],NULL);
    }
  }

  for (int t = NTHREADS-1; t >= 0; t--)
    Free_Profiles(paramc[t].P);
  free(paramc);
  for (int c = 0; c < N_OTYPE; c++)
    { for (int t = 0; t < NTHREADS; t++)
        free(paramm[c].fnames[t]);
      free(paramm[c].fnames);
      free(paramm[c].final);
    }
  free(paramm);
  free(path);
  free(root);
  free(FK_ROOT);
  free(threads);
  free_emodel(emodel);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  return 0;
}
