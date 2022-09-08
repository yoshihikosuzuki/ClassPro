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

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static char *Usage = "<relative_profile>[.prof] <source>[.db|.dam|.f[ast][aq][.gz] ";

#define N_EXT 10
char *EXT[N_EXT] = { ".db", ".dam",
                     ".fastq", ".fasta", ".fq", ".fa",
                     ".fastq.gz", ".fasta.gz", ".fq.gz", ".fa.gz" };

int main(int argc, char *argv[])
{ Profile_Index *P;

  // fastx
  gzFile      fxfp = NULL;
  kseq_t     *fxseq = NULL;

  // dazz_db
  DAZZ_DB    _db, *db = &_db;
  DAZZ_READ  *r = NULL;
  DAZZ_STUB  *stub = NULL;
  char      **flist = NULL;
  int        *findx = NULL;
  int         map;
  FILE       *hdrs = NULL;
  char       *hdrs_name = "";
  bool        is_db, is_dam;

  FILE       *cfile;

  // Parse arguments
  { int    i, j, k;
    int    flags[128];

    (void) flags;

    char  *path, *root;
    int    fid, idx;

    ARG_INIT("prof2class");

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

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    root  = Root(argv[1],".prof");
    path  = PathTo(argv[1]);
    cfile = Fopen(Catenate(path,"/",root,".class"),"w");

    path = PathTo(argv[2]);
    for (idx = 0; idx < N_EXT; idx++)
      { root  = Root(argv[2],EXT[idx]);
        fid   = open(Catenate(path,"/",root,EXT[idx]),O_RDONLY);
        if (fid >= 0) break;
        free(root);
      }
    if (idx == N_EXT || fid < 0)
      { fprintf(stderr,"Cannot open %s as a .db|.dam or .f{ast}[aq][.gz] file\n",argv[2]);
        exit(1);
      }
    close(fid);

    is_db  = (idx <= 1);
    is_dam = (idx == 1);

    free(root);
    free(path);
  }

  // Open files
  { if ((P = Open_Profiles(argv[1])) == NULL)
      { fprintf(stderr,"%s: Cannot open %s as a .prof file\n",Prog_Name,argv[1]);
        exit (1);
      }

    if (is_db)
      { char *pwd, *root;
        int status;

        status = Open_DB(argv[2],db);
        if (status < 0)
          { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[2]);
            exit (1);
          }
        if (db->part > 0)
          { fprintf(stderr,"%s: Cannot be called on a block\n",Prog_Name);
            exit (1);
          }
        if (P->nreads != db->nreads)
          { fprintf(stderr,"Inconsistent # of reads: .prof (%d) != .db (%d)\n",P->nreads,db->nreads);
            exit(1);
          }

        if (!is_dam)
          { root      = Root(argv[2],".db");
            pwd       = PathTo(argv[2]);
            stub      = Read_DB_Stub(Catenate(pwd,"/",root,".db"),DB_STUB_NREADS|DB_STUB_PROLOGS);
            flist     = stub->prolog;
            findx     = stub->nreads;
            findx[-1] = 0;
            map       = 0;
          }
        else
          { root      = Root(argv[2],".dam");
            pwd       = PathTo(argv[2]);
            hdrs_name = Strdup(Catenate(pwd,"/.",root,".hdr"),"Allocating header file name");
            hdrs      = Fopen(hdrs_name,"r");
            if (hdrs_name == NULL || hdrs == NULL)
              { fprintf(stderr,"Cannot open .hdr file %s [errno=%d]\n",hdrs_name,errno);
                exit (1);
              }
            free(hdrs_name);
          }

        free(root);
        free(pwd);
      }
    else
      { fxfp = gzopen(argv[2], "r");
        if (fxfp == NULL)
          { fprintf(stderr,"%s: Cannot open %s [errno=%d]\n",Prog_Name,argv[2],errno);
            exit (1);
          }
        fxseq = kseq_init(fxfp);
      }
  }

  // Main process
  { int     id;
    uint16 *profile;
    int     pmax, plen;
    int     rlen, rlen_max;
    char   *seq = NULL;
    char   *asgn;
    char    header[MAX_NAME];

    const int Km1 = P->kmer-1;

    if (is_db)
      { rlen_max = db->maxlen;
        seq      = New_Read_Buffer(db);
      }
    else
      { rlen_max = 60000;   // FIXME: auto computation
      }

    pmax     = 20000;
    profile  = Malloc(pmax*sizeof(uint16),"Profile array");

    asgn = Malloc((rlen_max+1)*sizeof(char),"asgn");
    for (int i = 0; i < Km1; i++)
      asgn[i] = 'N';

    for (id = 0; id < P->nreads; id++)
      { if (is_db)
          { r = db->reads+id;
            rlen = r->rlen;
            Load_Read(db,id,seq,2);
          }
        else
          { kseq_read(fxseq);
            rlen = fxseq->seq.l;
            seq = (char *)fxseq->seq.s;
          }

        const int slen = strlen(seq);
        if (rlen != slen)
          { fprintf(stderr,"rlen (%d) != strlen(seq) (%d)\n",rlen,slen);
            exit(1);
          }
        if (rlen > rlen_max)
          { fprintf(stderr,"rlen (%d) > rlen_max (%d)\n",rlen,rlen_max);
            exit(1);
          }

        plen = Fetch_Profile(P,(int64)id,pmax,profile);
        if (plen > pmax)
          { pmax    = 1.2*plen + 1000;
            profile = Realloc(profile,pmax*sizeof(uint16),"Profile array");
            Fetch_Profile(P,(int64)id,pmax,profile);
          }

        if (is_db)
          { if (!is_dam)
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
        else
          { sprintf(header,"@%s %s",fxseq->name.s,fxseq->comment.s);
          }

        if (rlen <= Km1)
          { asgn[rlen] = '\0';
            fprintf(cfile,"%s\n%s\n+\n%s\n",header,seq,asgn);
            asgn[rlen] = 'N';
            continue;
          }

        int idx = Km1;
        for (int i = 0; i < plen; i++)
          { char c;
            switch (profile[i])
              { case 0:
                  c = 'E';
                  break;
                case 1:
                  c = 'H';
                  break;
                case 2:
                  c = 'D';
                  break;
                default:
                  c = 'R';
                  break;
              }
            asgn[idx++] = c;
          }
        asgn[idx] = '\0';

        fprintf(cfile,"%s\n%s\n+\n%s\n",header,seq,asgn);
      }
    free(profile);
  }

  // Epilog
  { fclose(cfile);
    Free_Profiles(P);

    if (is_db)
      { Close_DB(db);
        if (!is_dam)
          Free_DB_Stub(stub);
        else
          fclose(hdrs);
      }
    else
      { kseq_destroy(fxseq);
        gzclose(fxfp);
      }

    Catenate(NULL,NULL,NULL,NULL);
    Numbered_Suffix(NULL,0,NULL);
    free(Prog_Name);
  }

  return 0;
}
