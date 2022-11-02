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

void *merge_anno(void *arg)
{ Merge_Arg  *data     = (Merge_Arg *)arg;
  char      **onames   = data->onames;
  char       *ofinal   = data->ofinal;
  int         N        = data->nfiles;

  FILE      *f, *g;
  int64      offset, idx;

  f = Fopen(onames[0],"ab+");
  if (f == NULL)
    { fprintf(stderr,"Cannot open %s [errno=%d]\n",onames[0],errno);
      exit(1);
    }

  if (fseek(f,-(long)sizeof(int64),SEEK_END) == -1)
    { fprintf(stderr,"Skip header failed for %s [errno=%d]\n",onames[0],errno);
      exit(1);
    }

  int ret = fread(&offset,sizeof(int64),1,f);
  if (ret != 1)
    { fprintf(stderr,"Cannot read last index of the first .anno file [ret=%d]\n",ret);
      exit(1);
    }

  for (int i = 1; i < N; i++)
    { g = Fopen(onames[i],"rb");
      if (g == NULL)
        { fprintf(stderr,"Cannot open %s [errno=%d]\n",onames[i],errno);
          exit(1);
        }

      while (fread(&idx,sizeof(int64),1,g) == 1)
        { idx += offset;
          fwrite(&idx,sizeof(int64),1,f);
        }
      offset = idx;

      fclose(g);
      unlink(onames[i]);
    }
  fclose(f);

  char *command = Malloc(sizeof(char)*(strlen(onames[0])+strlen(ofinal)+10),"Allocating command");
  sprintf(command,"mv %s %s",onames[0],ofinal);
  if (system(command) == -1)
    { fprintf(stderr,"Cannot mv %s to %s\n",onames[0],ofinal);
      exit(1);
    }
  free(command);

  return (NULL);
}

void *merge_files(void *arg)
{ Merge_Arg  *data      = (Merge_Arg *)arg;
  char      **onames    = data->onames;
  char       *ofinal    = data->ofinal;
  int         N         = data->nfiles;
  bool        is_bin    = data->is_bin;

  FILE *f, *g;
  char *buf = Malloc(sizeof(char)*MERGE_BUF_SIZE,"Allocating merge buffer");
  int n;

  f = Fopen(onames[0], is_bin ? "ab+" : "a+");
  if (f == NULL)
    { fprintf(stderr,"Cannot open %s [errno=%d]\n",onames[0],errno);
      exit(1);
    }

  for (int i = 1; i < N; i++)
    { g = Fopen(onames[i], is_bin ? "rb" : "r");
      if (g == NULL)
        { fprintf(stderr,"Cannot open %s [errno=%d]\n",onames[i],errno);
          exit(1);
        }

      while ((n = fread(buf,sizeof(char),MERGE_BUF_SIZE,g)) > 0)
        fwrite(buf,sizeof(char),n,f);

      fclose(g);
      unlink(onames[i]);
    }
  fclose(f);
  free(buf);

  char *command = Malloc(sizeof(char)*(strlen(onames[0])+strlen(ofinal)+10),"Allocating command");
  sprintf(command,"mv %s %s",onames[0],ofinal);
  if (system(command) == -1)
    { fprintf(stderr,"Cannot mv %s to %s\n",onames[0],ofinal);
      exit(1);
    }
  free(command);

  return (NULL);
}

static inline int ndigit(int n)
{ if (n == 0)
    return 1;
  int d = (n >= 0) ? 0 : 1;
  for (; n != 0; d++)
    n /= 10;
  return d;
}

static void prepare_db(Arg *arg, Class_Arg *paramc, Merge_Arg *paramm)
{ // NOTE: Assuming there is only a single input .db or .dam file
  char     *name   = arg->snames[0];
  char     *path   = PathTo(name);
  char     *root   = Root(name,NULL);
  
  // 14 for "/." & ".class.xxxx."; 10 for thread ID & '\0'
  const int MAX_FN = MAX(strlen(path),strlen(arg->tmp_path))+strlen(root)+14+10;
  
  // Set output file names (for both scattering and merging)
  for (int c = 0; c < N_OTYPE; c++)
    { Out_Info o = O_INFO[c];
      
      paramm[c].onames = Malloc(sizeof(char *)*arg->nthreads,"Allocating fnames");
      for (int t = 0; t < arg->nthreads; t++)
        { paramm[c].onames[t] = Malloc(sizeof(char)*MAX_FN,"Allocating fname");
          sprintf(paramm[c].onames[t],"%s%s%s%s.%d",arg->tmp_path,o.sep,arg->out_root,o.suf,t+1);
        }

      paramm[c].ofinal = Malloc(sizeof(char)*MAX_FN,"Allocating fname");
      sprintf(paramm[c].ofinal,"%s%s%s%s",path,o.sep,arg->out_root,o.suf);

      paramm[c].nfiles = arg->nthreads;
      paramm[c].is_bin = o.is_bin;
    }

  // Set variables used during multi-thread classification
  for (int t = 0; t < arg->nthreads; t++)
    { paramc[t].db = Malloc(sizeof(DAZZ_DB),"Allocating dazz db");
      if (Open_DB(name,paramc[t].db) < 0)
        { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,name);
          exit (1);
        }
      if (paramc[t].db->part > 0)
        { fprintf(stderr,"%s: Cannot be called on a block\n",Prog_Name);
          exit (1);
        }
      if (arg->nreads != paramc[t].db->nreads)
        { fprintf(stderr,"Inconsistent # of reads: .prof (%d) != .db (%d)\n",arg->nreads,paramc[t].db->nreads);
          exit(1);
        }
      if (!arg->is_dam)
        paramc[t].stub = Read_DB_Stub(name,DB_STUB_NREADS|DB_STUB_PROLOGS);
      else
        { char *hdrs_name = Strdup(Catenate(path,"/.",root,".hdr"),"Allocating header file name");
          paramc[t].hdrs  = Fopen(hdrs_name,"r");
          if (hdrs_name == NULL || paramc[t].hdrs == NULL)
            { fprintf(stderr,"Cannot open .hdr file %s [errno=%d]\n",hdrs_name,errno);
              exit (1);
            }
          free(hdrs_name);
        }
      paramc[t].beg = t*arg->nparts;
      paramc[t].end = MIN((t+1)*arg->nparts,arg->nreads);

#ifndef DEBUG_SINGLE
      paramc[t].afile = Fopen(paramm[ANNO].onames[t],"wb");
      paramc[t].dfile = Fopen(paramm[DATA].onames[t],"wb");
      if (paramc[t].afile == NULL || paramc[t].dfile == NULL)
        { fprintf(stderr,"Cannot open .*.class.*.%d\n",t+1);
          exit (1);
        }
      paramc[t].ranno = Fopen(paramm[RANNO].onames[t],"wb");
      paramc[t].rdata = Fopen(paramm[RDATA].onames[t],"wb");
      if (paramc[t].ranno == NULL || paramc[t].rdata == NULL)
        { fprintf(stderr,"Cannot open .*.rep.*.%d\n",t+1);
          exit (1);
        }
#endif
    }

#ifdef DEBUG_SINGLE
  return;
#endif

  // Set output file pointer to *.class file per thread
#ifndef PARALLEL_WRITE
  for (int t = 0; t < arg->nthreads; t++)
    { paramc[t].cfile = Fopen(paramm[CLASS].onames[t],"w");
      if (paramc[t].cfile == NULL)
        { fprintf(stderr,"Cannot open %s\n",paramm[CLASS].onames[t]);
          exit (1);
        }
    }
#else
  DAZZ_READ  *r;
  DAZZ_STUB  *stub = NULL;
  char      **flist = NULL;
  int        *findx = NULL;
  int         map;
  FILE       *hdrs = NULL;
  char       *hdrs_name = NULL;
  char        header[MAX_NAME];

  int64       csize;
  int64       coffset[arg->nthreads];
  int         id;

  if (!arg->is_dam)
    { stub      = Read_DB_Stub(name,DB_STUB_NREADS|DB_STUB_PROLOGS);
      flist     = stub->prolog;
      findx     = stub->nreads;
      findx[-1] = 0;
      map       = 0;
    }
  else
    { hdrs_name = Strdup(Catenate(path,"/.",root,".hdr"),"Allocating header file name");
      hdrs      = Fopen(hdrs_name,"r");
      if (hdrs_name == NULL || hdrs == NULL)
        { fprintf(stderr,"Cannot open .hdr file %s [errno=%d]\n",hdrs_name,errno);
          exit (1);
        }
    }

  // Total file size and write start position per thread
  csize = 0;
  id    = 0;
  for (int t = 0; t < arg->nthreads; t++)
    { coffset[t] = csize;
      for (int i = 0; i < arg->nparts && id < arg->nreads; i++, id++)
        { r = paramc[0].db->reads+id;
          if (!arg->is_dam)
            { while (id < findx[map-1])
                map -= 1;
              while (id >= findx[map])
                map += 1;
              sprintf(header,"@%s/%d/%d_%d",flist[map],r->origin,r->fpulse,r->fpulse+r->rlen);
            }
          else
            { // NOTE: '>' is included for .dam
              FSEEKO(hdrs,r->coff,SEEK_SET)
              FGETS(header,MAX_NAME,hdrs)
              header[strlen(header)-1] = '\0';
            }
          csize += 2*(r->rlen)+strlen(header)+6;   // NOTE: 6 bytes for {\n,\n,+,\n,\n}
        }
    }
  if (!arg->is_dam)
    Free_DB_Stub(stub);
  else
    { fclose(hdrs);
      free(hdrs_name);
    }

  // Allocate file size
  int fd = open(paramm[CLASS].ofinal,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
  if (fd == -1)
    { fprintf(stderr,"Cannot open/create %s [errono=%d]\n",paramm[CLASS].ofinal,errno);
      exit(1);
    }
  if (VERBOSE)
    { fprintf(stderr,"Allocating %lld bytes for %s\n",csize,paramm[CLASS].ofinal);
      fflush(stderr);
    }
  int ret = posix_fallocate(fd,0,sizeof(char)*csize);   // TODO: macOS
  if (ret != 0)
    { fprintf(stderr,"fallocate failed [ret=%d]\n",ret);
      exit(1);
    }
  close(fd);

  // Move file pointer to the offset for each thread
  for (int t = 0; t < arg->nthreads; t++)
    { paramc[t].cfile = Fopen(paramm[CLASS].ofinal,"w");
      if (paramc[t].cfile == NULL)
        { fprintf(stderr,"Cannot open %s [errono=%d]\n",paramm[CLASS].ofinal,errno);
          exit(1);
        }
      if (fseek(paramc[t].cfile,sizeof(char)*coffset[t],SEEK_SET) == -1)
        { fprintf(stderr,"fseek failed [errono=%d]\n",errno);
          exit(1);
        }
    }
#endif   // PARALLEL_WRITE

  // Write header info to .anno.1
  { const int64 idx = 0;
    int size;

    // Data track for classifications/seeds
    size = 8;
    fwrite(&arg->nreads,sizeof(int),1,paramc[0].afile);
    fwrite(&size,sizeof(int),1,paramc[0].afile);
    fwrite(&idx,sizeof(int64),1,paramc[0].afile);
    
    // Mask track for repeats
    size = 0;
    fwrite(&arg->nreads,sizeof(int),1,paramc[0].ranno);
    fwrite(&size,sizeof(int),1,paramc[0].ranno);
    fwrite(&idx,sizeof(int64),1,paramc[0].ranno);
  }

  return;
}

static void prepare_fx(Arg *arg, Class_Arg *paramc, Merge_Arg *paramm)
{ // NOTE: Assuming there is only a single input FASTX file   // TODO: remove the limitation
  char     *name   = arg->snames[0];
  char     *path   = PathTo(name);
  char     *root   = Root(name,NULL);
  
  // 14 for "/." & ".class.xxxx."; 10 for thread ID & '\0'
  const int MAX_FN = MAX(strlen(path),strlen(arg->tmp_path))+strlen(root)+14+10;
  
  // Set output file names (for both scattering and merging)
  for (int c = 0; c < N_OTYPE; c++)
    { Out_Info o = O_INFO[c];
      
      paramm[c].onames = Malloc(sizeof(char *)*arg->nthreads,"Allocating fnames");
      for (int t = 0; t < arg->nthreads; t++)
        { paramm[c].onames[t] = Malloc(sizeof(char)*MAX_FN,"Allocating fname");
          sprintf(paramm[c].onames[t],"%s%s%s%s.%d",arg->tmp_path,o.sep,arg->out_root,o.suf,t+1);
        }

      paramm[c].ofinal = Malloc(sizeof(char)*MAX_FN,"Allocating fname");
      sprintf(paramm[c].ofinal,"%s%s%s%s",path,o.sep,arg->out_root,o.suf);

      paramm[c].nfiles = arg->nthreads;
      paramm[c].is_bin = o.is_bin;
    }

  // Set variables used during multi-thread classification
  for (int t = 0; t < arg->nthreads; t++)
    { paramc[t].fxfp  = gzopen(name, "r");
      if (paramc[t].fxfp == NULL)
        { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,name);
          exit (1);
        }
      paramc[t].fxseq = kseq_init(paramc[t].fxfp);
      
      paramc[t].beg = t*arg->nparts;
      paramc[t].end = MIN((t+1)*arg->nparts,arg->nreads);
    }

#ifdef DEBUG_SINGLE
  return;
#endif

  // Set output file pointer to *.class file per thread
#ifndef PARALLEL_WRITE
  for (int t = 0; t < arg->nthreads; t++)
    { paramc[t].cfile = Fopen(paramm[CLASS].onames[t],"w");
      if (paramc[t].cfile == NULL)
        { fprintf(stderr,"Cannot open %s\n",paramm[CLASS].onames[t]);
          exit (1);
        }
    }
#else
  fprintf(stderr,"Parallel write for FASTX is currently unsupported.\n");
  exit(1);
#endif   // PARALLEL_WRITE

  return;
}

void prepare_param(Arg *arg, Class_Arg *paramc, Merge_Arg *paramm)
{ if (arg->is_db)
    prepare_db(arg,paramc,paramm);
  else
    prepare_fx(arg,paramc,paramm);

  return;
}

void free_param(Arg *arg, Class_Arg *paramc, Merge_Arg *paramm)
{ for (int t = arg->nthreads-1; t >= 0; t--)
    { Free_Profiles(paramc[t].P);
      if (arg->is_db)
        { Close_DB(paramc[t].db);
          if (!arg->is_dam)
            Free_DB_Stub(paramc[t].stub);
          else
            fclose(paramc[t].hdrs);
        }
      else
        { kseq_destroy(paramc[t].fxseq);
          gzclose(paramc[t].fxfp);
        }
    }
  free(paramc);

  for (int c = 0; c < N_OTYPE; c++)
    { for (int t = 0; t < arg->nthreads; t++)
        free(paramm[c].onames[t]);
      free(paramm[c].onames);
      free(paramm[c].ofinal);
    }
  free(paramm);

  return;
}