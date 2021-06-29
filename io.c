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

#define BUF_SIZE 4096

void *merge_anno(void *arg)
{ Merge_Arg  *data     = (Merge_Arg *)arg;
  char      **fnames   = data->fnames;
  char       *final    = data->final;
  int         N        = data->N;

  FILE      *f, *g;
  int64      offset, idx;

  f = Fopen(fnames[0],"ab+");
  if (f == NULL)
    { fprintf(stderr,"Cannot open %s [errno=%d]\n",fnames[0],errno);
      exit(1);
    }

  if (fseek(f,-(long)sizeof(int64),SEEK_END) == -1)
    { fprintf(stderr,"Skip header failed for %s [errno=%d]\n",fnames[0],errno);
      exit(1);
    }

  int ret = fread(&offset,sizeof(int64),1,f);
  if (ret != 1)
    { fprintf(stderr,"Cannot read last index of the first .anno file [ret=%d]\n",ret);
      exit(1);
    }

  for (int i = 1; i < N; i++)
    { g = Fopen(fnames[i],"rb");
      if (g == NULL)
        { fprintf(stderr,"Cannot open %s [errno=%d]\n",fnames[i],errno);
          exit(1);
        }

      while (fread(&idx,sizeof(int64),1,g) == 1)
        { 
#ifdef DEBUG
          if (idx == 0)
            { fprintf(stderr,"First index must not be 0 except the first .anno file\n");
              exit(1);
            }
#endif

          idx += offset;
          fwrite(&idx,sizeof(int64),1,f);
        }
      offset = idx;

      fclose(g);
      unlink(fnames[i]);
    }
  fclose(f);
  
  if (rename(fnames[0],final) == -1)
    { fprintf(stderr,"Cannot rename %s to %s [errno=%d]\n",fnames[0],final,errno);
      exit(1);
    }

  return (NULL);
}

void *merge_files(void *arg)
{ Merge_Arg  *data      = (Merge_Arg *)arg;
  char      **fnames    = data->fnames;
  char       *final     = data->final;
  int         N         = data->N;
  bool        is_binary = data->is_binary;

  FILE *f, *g;
  char buf[BUF_SIZE];
  int n;

  f = Fopen(fnames[0], is_binary ? "ab+" : "a+");
  if (f == NULL)
    { fprintf(stderr,"Cannot open %s [errno=%d]\n",fnames[0],errno);
      exit(1);
    }

  for (int i = 1; i < N; i++)
    { g = Fopen(fnames[i], is_binary ? "rb" : "r");
      if (g == NULL)
        { fprintf(stderr,"Cannot open %s [errno=%d]\n",fnames[i],errno);
          exit(1);
        }

      while ((n = fread(buf,sizeof(char),BUF_SIZE,g)) > 0)
        fwrite(buf,sizeof(char),n,f);

      fclose(g);
      unlink(fnames[i]);
    }
  fclose(f);
  
  if (rename(fnames[0],final) == -1)
    { fprintf(stderr,"Cannot rename %s to %s [errno=%d]\n",fnames[0],final,errno);
      exit(1);
    }

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

void prepare_db(char *path, char *root, Class_Arg *paramc)
{ char *name = Catenate(path,"/",root,".db");

  for (int t = 0; t < NTHREADS; t++)
    { paramc[t].db = Malloc(sizeof(DAZZ_DB),"Allocating dazz db");
      if (Open_DB(name,paramc[t].db) < 0)
        { fprintf(stderr,"%s: Cannot open %s.db\n",Prog_Name,name);
          exit (1);
        }
      if (paramc[t].db->part > 0)
        { fprintf(stderr,"%s: Cannot be called on a block\n",Prog_Name);
          exit (1);
        }
      if (NREADS != paramc[t].db->nreads)
        { fprintf(stderr,"Inconsistent # of reads: .prof (%d) != .db (%d)\n",NREADS,paramc[t].db->nreads);
          exit(1);
        }
      paramc[t].stub = Read_DB_Stub(name,DB_STUB_NREADS|DB_STUB_PROLOGS);
    }

#ifdef PARALLEL_WRITE
  DAZZ_STUB  *stub;
  DAZZ_READ  *r;
  char      **flist;
  int        *findx;
  int         map;
  char       *cfname;
  int64       csize;
  int64       coffset[NTHREADS];
  int         hsize, id;

  stub      = Read_DB_Stub(name,DB_STUB_NREADS|DB_STUB_PROLOGS);
  flist     = stub->prolog;
  findx     = stub->nreads;
  findx[-1] = 0;
  map       = 0;

  // Total file size and write start position per thread
  csize = 0;
  id    = 0;
  for (int t = 0; t < NTHREADS; t++)
    { coffset[t] = csize;
      for (int i = 0; i < NPARTS && id < NREADS; i++, id++)
        { while (id < findx[map-1])
            map -= 1;
          while (id >= findx[map])
            map += 1;

          r      = paramc[0].db->reads+id;
          hsize  = strlen(flist[map])+ndigit(r->origin)+ndigit(r->fpulse)+ndigit(r->fpulse+r->rlen);
          csize += 2*(r->rlen)+hsize+9;   // NOTE: 9 bytes for {@,/,/,_,\n,\n,+,\n,\n}
        }
    }
  Free_DB_Stub(stub);

  // Allocate file size
  cfname = Catenate(path,"/",root,".class");
  int fd = open(cfname,O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
  if (fd == -1)
    { fprintf(stderr,"open/create fail [errono=%d]\n",errno);
      exit(1);
    }
  if (VERBOSE)
    { fprintf(stderr,"Allocating %lld bytes for %s\n",csize,cfname);
      fflush(stderr);
    }
  int ret = posix_fallocate(fd,0,sizeof(char)*csize);
  if (ret != 0)
    { fprintf(stderr,"fallocate fail [ret=%d]\n",ret);
      exit(1);
    }
  close(fd);

  // Prepare file descripter + offset per thread
  for (int t = 0; t < NTHREADS; t++)
    { paramc[t].cfile = Fopen(cfname,"w");
      if (paramc[t].cfile == NULL)
        { fprintf(stderr,"open fail [errono=%d]\n",errno);
          exit(1);
        }
      if (fseek(paramc[t].cfile,sizeof(char)*coffset[t],SEEK_SET) == -1)
        { fprintf(stderr,"fseek fail [errono=%d]\n",errno);
          exit(1);
        }
    }

  free(cfname);
#endif

  return;
}

void prepare_fx(char *fnames[], int nfiles, Class_Arg *paramc, char *path, char *root)
{
  return;
}
