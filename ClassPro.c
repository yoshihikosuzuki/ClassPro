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

bool VERBOSE;
int  READ_LEN;
bool IS_DB;

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

static Arg *parse_arg(int argc, char *argv[])
{ Arg   *arg = Malloc(sizeof(Arg),"Allocating Arg");
  int    i, j, k;
  int    flags[128];
  char  *eptr;
  (void) flags;

  ARG_INIT("ClassPro");

  arg->nthreads = DEFAULT_NTHREADS;
  arg->cov      = -1;
  arg->rlen     = DEFAULT_RLEN;
  arg->tmp_path = DEFAULT_TMP_PATH;
  arg->fk_root  = NULL;
  
  j = 1;
  for (i = 1; i < argc; i++)
    if (argv[i][0] == '-')
      switch (argv[i][1])
      { default:
          ARG_FLAGS("vs")
          break;
        case 'T':
          ARG_POSITIVE(arg->nthreads,"Number of threads")
          break;
        case 'c':
          ARG_NON_NEGATIVE(arg->cov,"Estimated k-mer coverage")
          break;
        case 'r':
          ARG_POSITIVE(arg->rlen,"Average read length")
          break;
        case 'N':
          arg->fk_root = argv[i]+2;
          break;
        case 'P':
          arg->tmp_path = argv[i]+2;
          break;
      }
    else
      argv[j++] = argv[i];
  argc = j;

  if (argc < 2)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  arg->verbose    = flags['v'];
  arg->find_seeds = flags['s'];
  arg->nfiles     = argc-1;
  
  { int    fid, idx;
    char  *path, *root;

    path = PathTo(argv[1]);

    for (idx = 0; idx < N_EXT; idx++)
      { root  = Root(argv[1],EXT[idx]);
        fid   = open(Catenate(path,"/",root,EXT[idx]),O_RDONLY);
        if (fid >= 0) break;
        free(root);
      }
    if (idx == N_EXT || fid < 0)
      { fprintf(stderr,"Cannot open %s as a .db|.dam or .f{ast}[aq][.gz] file\n",arg->snames[0]);
        exit(1);
      }
    close(fid);

    arg->is_db = (idx <= 1) ? true : false;

    if (arg->is_db && arg->nfiles != 1)
      { fprintf(stderr,"Only single file is accepted for .db and .dam\n");
        exit(1);
      }

    if (arg->fk_root == NULL)
      arg->fk_root = Strdup(Catenate(path,"/",root,""),"Set fk_root");

    if (arg->verbose)
      fprintf(stderr,"%d seq files (first file: path = %s  root = %s  ext = %s)\nFASTK root = %s\n",
                     arg->nfiles,path,root,EXT[idx],arg->fk_root);

    // Set full path
    arg->snames = Malloc(sizeof(char *)*arg->nfiles,"Allocating fnames");
    for (int i = 0; i < arg->nfiles; i++)
      { path = PathTo(argv[i+1]);
        root = Root(argv[i+1],EXT[idx]);
        arg->snames[i] = Strdup(Catenate(path,"/",root,EXT[idx]),"Set fname");
        if ((fid = open(arg->snames[i],O_RDONLY)) == -1)
          { fprintf(stderr,"Cannot open %s [errno=%d]\n",arg->snames[i],errno);
            exit(1);
          }
        close(fid);
      }
  }

  { char  *cpath, *spath;
    DIR   *dirp;

    if (arg->tmp_path[0] != '/')
      { cpath = getcwd(NULL,0);
        if (arg->tmp_path[0] == '.')
          { if (arg->tmp_path[1] == '/')
              spath = Catenate(cpath,arg->tmp_path+1,"","");
            else if (arg->tmp_path[1] == '\0')
              spath = cpath;
            else
              { fprintf(stderr,"\n%s: -P option: . not followed by /\n",Prog_Name);
                exit (1);
              }
          }
        else
          spath = Catenate(cpath,"/",arg->tmp_path,"");
        arg->tmp_path = Strdup(spath,"Allocating path");
        free(cpath);
      }
    else
      arg->tmp_path = Strdup(arg->tmp_path,"Allocating path");

    if ((dirp = opendir(arg->tmp_path)) == NULL)
      { fprintf(stderr,"\n%s: -P option: cannot open directory %s\n",Prog_Name,arg->tmp_path);
        exit (1);
      }
    closedir(dirp);
    
    if (arg->verbose)
      fprintf(stderr,"Temp dir path = %s\n",arg->tmp_path);
  }

  return arg;
}

static void free_arg(Arg *arg)
{ for (int i = 0; i < arg->nfiles; i++)
    free(arg->snames[i]);
  free(arg->snames);
  free(arg->fk_root);
  free(arg->tmp_path);
  free(arg);

  return;
}

int main(int argc, char *argv[])
{ Arg           *arg     = parse_arg(argc,argv);;
  Class_Arg     *paramc  = Malloc(sizeof(Class_Arg)*arg->nthreads,"Allocating class args");;
  Merge_Arg     *paramm  = Malloc(sizeof(Merge_Arg)*N_OTYPE,"Allocate merge args");;
  THREAD        *threads = Malloc(sizeof(THREAD)*arg->nthreads,"Allocating class threads");
  
  VERBOSE  = arg->verbose;
  READ_LEN = arg->rlen;
  
  Profile_Index *P;
  Error_Model   *emodel;

  // Precompute etc.
  { if ((P = Open_Profiles(arg->fk_root)) == NULL)
      { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,arg->fk_root);
        exit (1);
      }
    arg->nreads = P->nreads;
    arg->nparts = (arg->nreads / arg->nthreads) + (arg->nreads % arg->nthreads == 0 ? 0 : 1);

    if (arg->verbose)
      fprintf(stderr,"Total # of reads = %d, # of reads per thread = %d\n",arg->nreads,arg->nparts);

    precompute_probs();
    process_global_hist(arg->fk_root,arg->cov);
    emodel = calc_init_thres();
  }

  // Set file names, file pointers, etc. for each thread
  { for (int t = 0; t < arg->nthreads; t++)
      { paramc[t].emodel = emodel;

#ifndef DUP_PROFILE
        paramc[t].P = (t == 0) ? P : Clone_Profiles(P);
#else
        paramc[t].P = (t == 0) ? P : Open_Profiles(arg->fk_root);
#endif
        if (paramc[t].P == NULL)
          { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,arg->fk_root);
            exit (1);
          }
      }

    prepare_param(arg,paramc,paramm);
  }

  // Classification
  { if (arg->verbose)
      fprintf(stderr,"Classifying %d-mers%s...\n",P->kmer,arg->find_seeds ? "" : " & Finding seeds");

    for (int t = 1; t < arg->nthreads; t++)
      pthread_create(threads+t,NULL,kmer_class_thread,paramc+t);
    kmer_class_thread(paramc);
    for (int t = 1; t < arg->nthreads; t++)
      pthread_join(threads[t],NULL);
  }

  // Merging intermediate files into final output
  { int i;
#ifndef PARALLEL_WRITE
    const int b = CLASS;
#else
    const int b = CLASS+1;
#endif
    const int e = (arg->is_db) ? N_OTYPE : CLASS+1;

    if (VERBOSE)
      fprintf(stderr,"\nMerging files...\n");

    for (i = 0; i+1 < arg->nthreads && b+i < e-1; i++)
      { if (oann[b+i])
          pthread_create(threads+i+1,NULL,merge_anno,paramm+b+i);
        else
          pthread_create(threads+i+1,NULL,merge_files,paramm+b+i);
      }
    for (; b+i < e; i++)
      { if (oann[b+i])
          merge_anno(paramm+b+i);
        else
          merge_files(paramm+b+i);
      }
    for (i = 0; i+1 < arg->nthreads && b+i < e-1; i++)
      pthread_join(threads[i+1],NULL);
  }

  // Epilogue
  { free(threads);
    free_emodel(emodel);
    free_param(arg,paramc,paramm);
    free_arg(arg);

    Catenate(NULL,NULL,NULL,NULL);
    Numbered_Suffix(NULL,0,NULL);
    free(Prog_Name);
  }

  return 0;
}
