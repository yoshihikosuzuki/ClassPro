#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "libfastk.h"
#include "DB.h"

static char *Usage = "<source_root> <E/H_thres> <H/D_thres> <D/R_thres>";

int main(int argc, char *argv[])
{ Profile_Index *P;

  // NOTE: Only DAZZ_DB is accepted for now
  DAZZ_DB    _db, *db = &_db;
  DAZZ_STUB  *stub = NULL;
  char      **flist = NULL;
  int        *findx = NULL;
  int         map;
  FILE       *hdrs = NULL;
  char       *hdrs_name = "";
  bool        is_dam;

  FILE       *cfile;

  int  THRES[4];

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("ClassGS");

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

    if (argc != 5)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    THRES[0]= (int)strtol(argv[2],&eptr,10);
    THRES[1]= (int)strtol(argv[3],&eptr,10);
    THRES[2]= (int)strtol(argv[4],&eptr,10);

    fprintf(stderr,"E < %d <= H < %d <= D < %d <= R\n",THRES[0],THRES[1],THRES[2]);
  }

  { char *pwd, *root;
    int status;

    status = Open_DB(argv[1],db);
    if (status < 0)
      { fprintf(stderr,"%s: Cannot open %s as .db/.dam\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (db->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block\n",Prog_Name);
        exit (1);
      }
    is_dam = (status == 1) ? true : false;

    if (!is_dam)
      { root      = Root(argv[1],".db");
        pwd       = PathTo(argv[1]);
        stub      = Read_DB_Stub(Catenate(pwd,"/",root,".db"),DB_STUB_NREADS|DB_STUB_PROLOGS);
        flist     = stub->prolog;
        findx     = stub->nreads;
        findx[-1] = 0;
        map       = 0;
      }
    else
      { root      = Root(argv[1],".dam");
        pwd       = PathTo(argv[1]);
        hdrs_name = Strdup(Catenate(pwd,"/.",root,".hdr"),"Allocating header file name");
        hdrs      = Fopen(hdrs_name,"r");
        if (hdrs_name == NULL || hdrs == NULL)
          { fprintf(stderr,"Cannot open .hdr file %s [errno=%d]\n",hdrs_name,errno);
            exit (1);
          }
        free(hdrs_name);
      }
    
    cfile = Fopen(Catenate(pwd,"/",root,".GS.class"),"w");

    free(root);
    free(pwd);
  }

  { P = Open_Profiles(argv[1]);
    if (P == NULL)
      { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (db->nreads != P->nreads)
      { fprintf(stderr,"# of reads must be same\n");
        exit(1);
      }
  }

  { int     id;
    DAZZ_READ *r;
    uint16 *profile;
    int     pmax, plen;
    int     rlen, rlen_max;
    char   *seq;
    char   *asgn;
    char    header[MAX_NAME];
    const int Km1 = P->kmer-1;

    rlen_max = db->maxlen;
    seq      = New_Read_Buffer(db);

    pmax     = 20000;
    profile  = Malloc(pmax*sizeof(uint16),"Profile array");

    asgn = Malloc((rlen_max+1)*sizeof(char),"buf");
    for (int i = 0; i < Km1; i++)
      asgn[i] = 'N';

    for (id = 0; id < P->nreads; id++)
      { r = db->reads+id;
        rlen = r->rlen;
        Load_Read(db,id,seq,2);

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

        if (!is_dam)
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
        
        if (rlen <= Km1)
          { asgn[rlen] = '\0';
            fprintf(cfile,"%s\n%s\n+\n%s\n",header,seq,asgn);
            asgn[rlen] = 'N';
            continue;
          }

        int idx = Km1;
        for (int i = 0; i < plen; i++)
          { if (profile[i] < THRES[0])
              asgn[idx++] = 'E';
            else if (profile[i] < THRES[1])
              asgn[idx++] = 'H';
            else if (profile[i] < THRES[2])
              asgn[idx++] = 'D';
            else
              asgn[idx++] = 'R';
          }
        asgn[idx] = '\0';

        fprintf(cfile,"%s\n%s\n+\n%s\n",header,seq,asgn);
      }
    free(profile);
  }

  fclose(cfile);

  Free_Profiles(P);
  Close_DB(db);
  if (!is_dam)
    Free_DB_Stub(stub);
  else
    fclose(hdrs);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}