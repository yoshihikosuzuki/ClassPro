#include <stdio.h>
#include <stdbool.h>
#include "libfastk.h"
#include "DB.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define WSIZE 1000

const char *Usage = "[-v] <fastk-prefix> <dazz_db>";

int main(int argc, char *argv[])
{ char *fk_root, *db_fname;
  bool  verbose;
  gzFile   F;
  kseq_t   *ks;
  DAZZ_DB _db, *db = &_db;
  Profile_Index *P;
  FILE *kanno, *kdata;
  FILE *ranno, *rdata;
  FILE *danno, *ddata;
  FILE *sanno, *sdata;
  int64 ridx, sidx;

  // Parse arguments
  { int    i, j, k;
    int    flags[128];
    // char  *eptr;
    (void) flags;

    ARG_INIT("RepDepAnno");

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    verbose = flags['v'];
    fk_root = argv[1];
    db_fname = argv[2];
  }

  // Load profiles and classifications
  { if (verbose)
      fprintf(stderr,"Opening %s.prof\n",fk_root);

    if ((P = Open_Profiles(fk_root)) == NULL)
      { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,fk_root);
        exit (1);
      }
    
    if (verbose)
      fprintf(stderr,"Opening %s\n",db_fname);
    int status = Open_DB(db_fname,db);
    if (status < 0)
      { fprintf(stderr,"%s: Cannot open %s as .db/.dam\n",Prog_Name,db_fname);
        exit (1);
      }
    Load_All_Reads(db,0);

    if (verbose)
      fprintf(stderr,"Opening %s.class\n",fk_root);

    if ((F = gzopen(Catenate(fk_root,"","",".class"), "r")) == NULL)
      { fprintf(stderr,"%s: Cannot open %s.class [errno=%d]\n",Prog_Name,fk_root,errno);
        exit (1);
      }
    ks = kseq_init(F);

    if (verbose)
      fprintf(stderr,"Opening .%s.class.[anno|data]\n",fk_root);

    kanno = Fopen(Catenate(".",fk_root,"",".class.anno"),"r");
    kdata = Fopen(Catenate(".",fk_root,"",".class.data"),"r");
    if (kanno == NULL || kdata == NULL)
      exit(1);
  }

  // Open output tracks
  { int size = 0;
    int64 idx = 0;
    
    // ranno = Fopen(Catenate(".",fk_root,"",".rep.anno"),"w");
    // rdata = Fopen(Catenate(".",fk_root,"",".rep.data"),"w");
    // if (ranno == NULL || rdata == NULL)
    //   exit(1);
    // fwrite(&(db->nreads),sizeof(int),1,ranno);
    // fwrite(&size,sizeof(int),1,ranno);
    // fwrite(&idx,sizeof(int64),1,ranno);
    // ridx = 0;

    sanno = Fopen(Catenate(".",fk_root,"",".seed.anno"),"w");
    sdata = Fopen(Catenate(".",fk_root,"",".seed.data"),"w");
    if (sanno == NULL || sdata == NULL)
      exit(1);
    fwrite(&(db->nreads),sizeof(int),1,sanno);
    fwrite(&size,sizeof(int),1,sanno);
    fwrite(&idx,sizeof(int64),1,sanno);
    sidx = 0;
  }

  // Main routine
  { uint16 *profile = Malloc(db->maxlen*sizeof(uint16),"Profile array");
    char *cls = New_Read_Buffer(db), *c;
    char *s = (char *) (db->bases);
    int rlen;
    uint32 clen;
    int64  o, p;

    o = 0;
    for (int64 id = 0; id < db->nreads; id++)
      { rlen  = db->reads[id].rlen;
        
        if (verbose)
          fprintf(stderr,"Read %lld (%5d bp):\n",id+1,rlen);
        
        if (verbose)
          { fprintf(stderr,"  POS = ");
            for (int i = 0; i < rlen; i++)
              fprintf(stderr," %5d",i);
            fprintf(stderr,"\n");
          }

        // Prof
        int plen = Fetch_Profile(P,id,db->maxlen,profile);
        if (verbose)
          { fprintf(stderr," PROF = ");
            for (int i = 0; i < plen; i++)
              fprintf(stderr," %5d",profile[i]);
            fprintf(stderr,"\n");
          }

        // Class
        if (kseq_read(ks) < 0)
          { fprintf(stderr,"Cannot load seq from %s.class\n",fk_root);
            exit(1);
          }
        { int Km1 = P->kmer-1;
          int b = Km1, e;
          bool in_R = (ks->qual.s[Km1] == 'R') ? true : false;
          int nrep = 0;
          for (int i = Km1+1; i < rlen; i++)
            { if (!in_R)
                { if (ks->qual.s[i] == 'R')
                    { b = i;
                      in_R = true;
                    }
                }
              if (in_R)
                { if (ks->qual.s[i] != 'R')
                    { e = i;
                      // fwrite(&b,sizeof(int),1,rdata);
                      // fwrite(&e,sizeof(int),1,rdata);
                      nrep++;
                      in_R = false;
                    }
                }
            }
          if (in_R)
            { e = rlen;
              // fwrite(&b,sizeof(int),1,rdata);
              // fwrite(&e,sizeof(int),1,rdata);
              nrep++;
            }
          ridx += (nrep*2*sizeof(int));
          // fwrite(&ridx,sizeof(int64),1,ranno);
        }
        if (verbose)
          { fprintf(stderr,"CLASS = ");
            for (int i = 0; i < rlen; i++)
              fprintf(stderr," %5c",ks->qual.s[i]);
            fprintf(stderr,"\n");
          }

        // Seed
        // TODO: Load class without DB? (know clen by anno)
        clen = COMPRESSED_LEN(rlen); 
        if (clen > 0)
          { if (fread(cls,1,clen,kdata) != clen)
              { fprintf(stderr,"%s: Cannot load read %lld\n",Prog_Name,id);
                exit (1);
              }
          }
        Uncompress_Read(rlen,cls);
        c = cls-o;
        int nseed = 0;
        for (p = o; s[p] != 4; p++)   // TODO: loop without DB? (length must be rlen)
          { if (c[p] > 0)
              { // int b = MAX(p-o-WSIZE,0);
                // int e = MIN(p-o+WSIZE,rlen);
                int b = p-o-P->kmer+1;
                if (b < 0)
                  { fprintf(stderr,"[ERROR] position < Km1\n");
                    exit(1);
                  }
                int e = p-o;
                fwrite(&b,sizeof(int),1,sdata);
                fwrite(&e,sizeof(int),1,sdata);
                nseed++;
              }
          }
        sidx += (nseed*2*sizeof(int));
        fwrite(&sidx,sizeof(int64),1,sanno);
        if (verbose)
          { fprintf(stderr,"SEED:\n");
            for (p = o; s[p] != 4; p++)
              if (c[p] > 0)
                fprintf(stderr,"idx = %lld, c[idx] = %d\n",p-o,c[p]);
          }

        o += (rlen+1);

        if (verbose)
          fprintf(stderr,"\n");
      }
  }

  fclose(kdata);
  // fclose(ranno);
  // fclose(rdata);
  fclose(sanno);
  fclose(sdata);

  return 0;
}