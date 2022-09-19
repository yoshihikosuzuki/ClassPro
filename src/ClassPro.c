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

#include "ClassPro.h"
#include "benchmark.h"

#include "const.c"
#include "io.c"
#include "prob.c"
#include "util.c"
#include "hist.c"
#include "context.c"
#include "wall.c"
#include "class_rel.c"
#include "class_unrel.c"
#include "seed.c"

bool  VERBOSE;
int   READ_LEN;
bool  IS_DB;
bool  IS_DAM;
bool  FIND_SEED;
cnt_t GLOBAL_COV[N_STATE];

static void *kmer_class_thread(void *arg)
{ // Shared parameters & constants
  Class_Arg     *data   = (Class_Arg *)arg;
  Profile_Index *P      = data->P;
  Error_Model   *emodel = data->emodel;
  const int      K      = P->kmer;
  const int      Km1    = K-1;
  kdq_t(hmer_t) *Q      = kdq_init(hmer_t);

#ifdef DEBUG_SINGLE
  if (DEBUG_SINGLE_ID < data->beg+1 || data->end < DEBUG_SINGLE_ID)
    goto class_exit;
#endif

  // Variables for classification   // TODO: reduce memory by reusing variables
  int       rlen_max;
  int       rlen, plen;        // `rlen` = `plen` + `Km1`
  char     *seq = NULL;        // Fetched sequence of the read
  cnt_t    *profile;           // Fetched profile of a read
  Seq_Ctx  *ctx[N_WTYPE];      // Context lengths per position
  Seq_Ctx  *lctx, *_lctx;
  Seq_Ctx  *rctx;
  int       N, M;              // Number of walls/reliable intervals
  Intvl    *intvl, *rintvl;
  char     *rasgn, *pasgn;     // Classifications of intervals or k-mers for read and profile
  int      *sasgn = NULL;
  int      *hash = NULL;
  seg_t   *cprofile = NULL;   // Compressed profile, used in seed selection
  intvl_t *mintvl = NULL;
  Wall_Arg *warg;
  Rel_Arg  *rel_arg;
#ifdef DO_PMM
  PMM_Arg  *parg;
#endif

  // Variables for loading sequence/header
  kseq_t     *fxseq = NULL;       // for FASTX input
  DAZZ_DB    *db = NULL;
  DAZZ_READ  *r = NULL;
#ifdef WRITE_TRACK
  int64       tidx;               // Index for DAZZ track (afile; .anno file)
  char       *track = NULL;       // Data for DAZZ track (dfile; .data file)
  char       *crack = NULL;
#endif
  char        header[MAX_NAME];   // for both FASTX and DB
  DAZZ_STUB  *stub;
  char      **flist = NULL;
  int        *findx = NULL;
  int         map;
  FILE       *hdrs = NULL;
  char       *hdrs_name = "";

  // Input type-specific processing
  if (IS_DB)
    { db       = data->db;
      rlen_max = db->maxlen;
      seq      = New_Read_Buffer(db);
#ifdef WRITE_TRACK
      tidx     = 0;
      track    = New_Read_Buffer(db);
      crack    = track + Km1;
#endif
      if (!IS_DAM)
        { stub  = data->stub;
          flist = stub->prolog;
          findx = stub->nreads;
          map   = 0;
          findx[-1] = 0;
        }
      else
        hdrs = data->hdrs;
    }
  else   // FASTX
    { fxseq = data->fxseq;
      for (int i = 0; i < data->beg; i++)
        if (kseq_read(fxseq) < 0)
          { fprintf(stderr,"Cannot load %d-th read\n",i+1);
            exit(1);
          }
      rlen_max = MAX_READ_LEN;
    }

  // Allocation
  { rasgn = Malloc((rlen_max+1)*sizeof(char),"Classification array");
    pasgn = rasgn + Km1;
    for (int i = 0; i < Km1; i++)
      rasgn[i] = 'N';
    if (FIND_SEED)
      { sasgn    = Malloc((rlen_max+1)*sizeof(int),"Seed array");
        hash     = Malloc((rlen_max+1)*sizeof(int),"Hash array");
        cprofile = Malloc((rlen_max+1)*sizeof(seg_t),"Compressed profile array");
        mintvl   = Malloc((rlen_max+1)*sizeof(intvl_t),"Mask interval array");
      }

    rel_arg = alloc_rel_arg(rlen_max);
    warg = alloc_wall_arg(rlen_max);
#ifdef DO_PMM
    parg = alloc_pmm_arg(rlen_max);
#endif

    intvl   = Malloc(rlen_max*sizeof(Intvl),"Interval array");
    rintvl  = Malloc(rlen_max*sizeof(Intvl),"Interval array");
    profile = Malloc(rlen_max*sizeof(cnt_t),"Profile array");

    _lctx   = Malloc(rlen_max*sizeof(Seq_Ctx),"Allocating left ctx vector");
    rctx    = Malloc(rlen_max*sizeof(Seq_Ctx),"Allocating right ctx vector");
    lctx    = _lctx + Km1 - 1;
    _lctx[0][HP] = 1;
    _lctx[0][DS] = _lctx[0][TS] = _lctx[1][TS] = 0;
    ctx[DROP] = lctx;
    ctx[GAIN] = rctx;
  }

  // For each read, do classification
  for (int id = data->beg; id < data->end; id++)
    {
#ifdef DEBUG_SINGLE
      if (id+1 < DEBUG_SINGLE_ID)
        { if (!IS_DB)
            kseq_read(fxseq);
          continue;
        }
      else if (DEBUG_SINGLE_ID < id+1)
        goto class_exit;
#endif

#ifdef WRITE_TRACK
      if (IS_DB)
        for (int i = 0; i < Km1; i++)
          track[i] = 0;
#endif

      // 1. Load sequence and header
      { if (IS_DB)
          { r = db->reads+id;
            rlen = r->rlen;
            Load_Read(db,id,seq,2);

            if (!IS_DAM)
              { while (id < findx[map-1])
                  map -= 1;
                while (id >= findx[map])
                  map += 1;
                sprintf(header,"@%s/%d/%d_%d",flist[map],r->origin,r->fpulse,r->fpulse+r->rlen);
              }
            else
              { FSEEKO(hdrs,r->coff,SEEK_SET)
                FGETS(header,MAX_NAME,hdrs)
                header[strlen(header)-1] = '\0';
                header[0] = '@';
              }
          }
        else   // FASTX
          { kseq_read(fxseq);
            rlen = fxseq->seq.l;
            seq = (char *)fxseq->seq.s;
            if (rlen > MAX_READ_LEN)   // TODO: realloc when exceeded
              { fprintf(stderr,"rlen (%d) > MAX_READ_LEN for FASTX inputs (%d)\n",rlen,MAX_READ_LEN);
                exit(1);
              }
            sprintf(header,"@%s %s",fxseq->name.s,fxseq->comment.s);
          }

#if defined(DEBUG_ITER) || defined(DEBUG_CTX) || defined(DEBUG_ERROR)
        fprintf(stderr,"\nRead %5d (%5d bp): ",id+1,rlen);
        fflush(stderr);
#endif

#ifdef DEBUG
        const int slen = strlen(seq);
        if (rlen != slen)
          { fprintf(stderr,"rlen (%d) != strlen(seq) (%d)\n",rlen,slen);
            exit(1);
          }
        if (rlen > rlen_max)
          { fprintf(stderr,"rlen (%d) > rlen_max (%d)\n",rlen,rlen_max);
            exit(1);
          }
#endif

        // Special case where the read is too short
        if (rlen <= Km1)
          {
#ifdef DEBUG_SINGLE
            fprintf(stderr,"rlen (%d) <= K-1 (%d)\n",rlen,Km1);
            continue;
#endif

            fprintf(data->cfile,"%s\n%s\n+\n%*s\n",header,seq,rlen,rasgn);

#ifdef WRITE_TRACK
            if (IS_DB)
              { Compress_Read(rlen,track);
                int t = COMPRESSED_LEN(rlen);
                fwrite(track,1,t,data->dfile);
                tidx += t;
                fwrite(&tidx,sizeof(int64),1,data->afile);
              }
#endif

            continue;
          }
      }

      // 2. Compute per-base sequence context
      calc_seq_context(_lctx,rctx,seq,rlen);

      // 3. Load count profile
      plen = Fetch_Profile(P,(int64)id,rlen_max,profile);
      if (rlen != plen+Km1)
        { fprintf(stderr,"Read %d: rlen (%d) != plen+Km1 (%d)\n",id+1,rlen,plen+Km1);
          exit(1);
        }

      // 4. Wall detection
      N = find_wall(warg,intvl,profile,plen,ctx,emodel,K);
#ifdef DEBUG_ITER
      fprintf(stderr,"%3d intvls",N);
#endif

      // 5. Reliable interval detection
      M = find_rel_intvl(intvl,N,rintvl,profile,ctx,K);
#ifdef DEBUG_ITER
      int rtot = 0;
      for (int i = 0; i < M; i++)
        rtot += rintvl[i].e - rintvl[i].b;
      fprintf(stderr,", %3d rel intvls (%%base=%2d)",M,(int)(100.*rtot/plen));
#endif

#ifdef DO_PMM
      // Local coverage estimation via Poisson mixture model (optional)
      int nnorm = pmm_vi(parg,profile,plen,lambda);
#ifdef DEBUG_ITER
      fprintf(stderr,", (H,D)=(%.lf,%.lf) (%2d %% normal)",lambda[0],lambda[1],(int)(100.*nnorm/plen));
#endif
#endif   // DO_PMM

      // 6. Classification of reliable intervals and then the rest
      { classify_rel(rel_arg,rintvl,M,intvl,N,plen);
        classify_unrel(intvl,N);
        for (int i = 0; i < N; i++)
          { Intvl I = intvl[i];
            char c = stoc[(int)I.asgn];
            for (pos_t j = I.b; j < I.e; j++)
              pasgn[j] = c;
          }
        rasgn[rlen] = '\0';
        // remove_slip(pasgn,profile,plen,ctx);

        // printf("r");
        // for (int i = 0; i < rlen; i++)
        //   printf("%c",rasgn[i]);
        // printf("\n");
      }

      // Find seeds for alignment
      if (FIND_SEED)
        find_seeds(Q,seq,pasgn,profile,cprofile,hash,sasgn,mintvl,plen,K);

#ifdef DEBUG_SINGLE
      continue;
#endif

      // 7. Output to .class file (and optionally DAZZ track)
      { fprintf(data->cfile,"%s\n%s\n+\n%s\n",header,seq,rasgn);
#ifdef WRITE_TRACK
        if (IS_DB)
          { for (int i = 0; i < plen; i++)
              crack[i] = (FIND_SEED) ? ctos[sasgn[i]] : ctos[(int)pasgn[i]];

            // printf("t");
            // for (int i = 0; i < rlen; i++)
            //   printf("%d",track[i]);
            // printf("\n");

            Compress_Read(rlen,track);
            int t = COMPRESSED_LEN(rlen);
            fwrite(track,1,t,data->dfile);
            tidx += t;
            fwrite(&tidx,sizeof(int64),1,data->afile);
          }
#endif
      }
    }   // Loop for each read

  fclose(data->cfile);
  if (IS_DB)
    { free(seq-1);
#ifdef WRITE_TRACK
      fclose(data->afile);
      fclose(data->dfile);
      free(track-1);
#endif
    }
  free(profile);
  free(_lctx);
  free(rctx);
  free(intvl);
  free(rintvl);
  free(rasgn);
  free_wall_arg(warg);
  free_rel_arg(rel_arg, rlen_max);
  kdq_destroy(hmer_t,Q);
#ifdef DO_PMM
  free_pmm_arg(parg);
#endif

#ifdef DEBUG_SINGLE
class_exit:
#endif

  return (NULL);
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

static Arg *parse_arg(int argc, char *argv[])
{ Arg   *arg = Malloc(sizeof(Arg),"Allocating Arg");
  int    i, j, k;
  int    flags[128];
  char  *eptr;
  (void) flags;

  ARG_INIT("ClassPro");

  arg->nthreads = (int)DEFAULT_NTHREADS;
  arg->cov      = 0;
  arg->rlen     = (int)DEFAULT_RLEN;
  arg->tmp_path = (char *)DEFAULT_TMP_PATH;
  arg->fk_root  = NULL;
  arg->out_root = NULL;
  arg->model_path  = NULL;

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
        case 'M':
          arg->model_path = argv[i]+2;
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
  FIND_SEED = arg->find_seeds = flags['s'];
  arg->nfiles     = argc-1;

  if (arg->verbose)
    fprintf(stderr,"Info about inputs:\n");

  // Input file names
  { int    fid, idx;
    char  *path, *root;

    path = PathTo(argv[1]);

    for (idx = 0; idx < N_EXT; idx++)
      { root  = Root(argv[1],EXT[idx]);
        fid   = open(Catenate(path,"/",root,EXT[idx]),O_RDONLY);
        if (fid >= 0)
          { arg->out_root = Strdup(root,"Set out_root");
            break;
          }
        free(root);
      }
    if (idx == N_EXT || fid < 0)
      { fprintf(stderr,"Cannot open %s as a .db|.dam or .f{ast}[aq][.gz] file\n",argv[1]);
        exit(1);
      }
    close(fid);

    IS_DB  = arg->is_db  = (idx <= 1) ? true : false;
    IS_DAM = arg->is_dam = (idx == 1) ? true : false;

    if (arg->is_db && arg->nfiles != 1)
      { fprintf(stderr,"Only single file is accepted for .db and .dam\n");
        exit(1);
      }

    // TODO: extend to multiple FASTX files
    if (arg->nfiles != 1)
      { fprintf(stderr,"Currently only single file is accepted for FASTX input\n");
        exit(1);
      }

    if (arg->fk_root == NULL)
      arg->fk_root = Strdup(Catenate(path,"/",root,""),"Set fk_root");

    if (arg->verbose)
      { fprintf(stderr,"    # of sequence files   = %d\n",arg->nfiles);
        fprintf(stderr,"    First (path,root,ext) = (%s, %s, %s)\n",path,root,EXT[idx]);
        fprintf(stderr,"    FASTK outputs' root   = %s\n",arg->fk_root);
        fprintf(stderr,"    Otput .class file     = %s/%s.class\n",path,arg->out_root);
      }

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

  // Temparary directory
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
      fprintf(stderr,"    Temp dir path         = %s\n",arg->tmp_path);
  }

  return arg;
}

int main(int argc, char *argv[])
{ startTime();

  Arg           *arg;
  Class_Arg     *paramc;
  Merge_Arg     *paramm;
  pthread_t     *threads;
  Profile_Index *P;
  Error_Model   *emodel;

  // Parse command-line arguments and allocate parameter objects
  { arg      = parse_arg(argc,argv);
    VERBOSE  = arg->verbose;
    READ_LEN = arg->rlen;

    paramc   = Malloc(sizeof(Class_Arg)*arg->nthreads,"Allocating class args");
    paramm   = Malloc(sizeof(Merge_Arg)*N_OTYPE,"Allocate merge args");
    threads  = Malloc(sizeof(pthread_t)*arg->nthreads,"Allocating class threads");
  }

  // Precompute several values
  { // 1. Number of reads processed per thread
    if ((P = Open_Profiles(arg->fk_root)) == NULL)
      { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,arg->fk_root);
        exit (1);
      }
    arg->nreads = P->nreads;
    arg->nparts = (arg->nreads / arg->nthreads) + (arg->nreads % arg->nthreads == 0 ? 0 : 1);
    if (arg->verbose)
      { fprintf(stderr,"    Total # of reads      = %d\n",arg->nreads);
        fprintf(stderr,"    # of reads per thread = %d\n",arg->nparts);
      }

    // 2. Constants for probability calculation
    precompute_logfact();

    // 3. Constants for Poisson mixture model
    // precompute_digamma();

    // 4. Global Haploid/Diploid coverages
    process_global_hist(arg->fk_root,arg->cov);
    GLOBAL_COV[HAPLO] = lambda_prior[0];
    GLOBAL_COV[DIPLO] = lambda_prior[1];
    GLOBAL_COV[ERROR] = 1;
    GLOBAL_COV[REPEAT] = plus_sigma(GLOBAL_COV[DIPLO],N_SIGMA_RCOV);
    DR_RATIO = 1.+(double)N_SIGMA_R*(1./sqrt(GLOBAL_COV[DIPLO]));
    if (arg->verbose)
      fprintf(stderr,"    Estimated R-threshold = %d\n",GLOBAL_COV[REPEAT]);

    // 5. Thresholds of a count change due to errors in self and errors in others
    //    Context-specific sequcning error model is also loaded
    emodel = calc_init_thres(arg->model_path);
  }

  // Set error model, profile reader, and output file names/pointers for each thread
  { for (int t = 0; t < arg->nthreads; t++)
      { paramc[t].emodel = emodel;
#ifdef DUP_PROFILE
        paramc[t].P = (t == 0) ? P : Open_Profiles(arg->fk_root);
#else
        paramc[t].P = (t == 0) ? P : Clone_Profiles(P);
#endif
        if (paramc[t].P == NULL)
          { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,arg->fk_root);
            exit (1);
          }
      }
    prepare_param(arg,paramc,paramm);
  }

  // Main classification routine
  { if (arg->verbose)
      fprintf(stderr,"Classifying %d-mers%s...\n",P->kmer,arg->find_seeds ? " & Finding seeds" : "");

    for (int t = 1; t < arg->nthreads; t++)
      pthread_create(threads+t,NULL,kmer_class_thread,paramc+t);
    kmer_class_thread(paramc);
    for (int t = 1; t < arg->nthreads; t++)
      pthread_join(threads[t],NULL);
  }

  if (arg->verbose)
    timeTo(stderr,false);

#ifdef DEBUG_SINGLE
  return 0;
#endif

  // Merging intermediate files into final output .class file
  { int i;
#ifndef PARALLEL_WRITE
    const int b = CLASS;
#else
    const int b = CLASS+1;   // Skip .class file because already merged
#endif
#ifdef WRITE_TRACK
    const int e = (arg->is_db) ? N_OTYPE : CLASS+1;
#else
    const int e = CLASS+1;
#endif

    if (VERBOSE && b < e)
      fprintf(stderr,"\nMerging files...\n");

    for (i = 0; i+1 < arg->nthreads && b+i < e-1; i++)
      { if (O_INFO[b+i].is_anno)
          pthread_create(threads+i+1,NULL,merge_anno,paramm+b+i);
        else
          pthread_create(threads+i+1,NULL,merge_files,paramm+b+i);
      }
    for (; b+i < e; i++)
      { if (O_INFO[b+i].is_anno)
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

  if (arg->verbose)
    timeTo(stderr,true);

  return 0;
}
