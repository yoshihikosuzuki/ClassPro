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

const char btoc[4] = { 'a', 'c', 'g', 't' };
const char ctos[4] = { 'E', 'R', 'H', 'D' };

const char *Usage = "[-v] <fastk-prefix> <dazz_db>";

int main(int argc, char *argv[])
{ char *fk_root, *db_fname;
  bool  verbose;
  gzFile   F;
  kseq_t   *ks;
  DAZZ_DB _db, *db = &_db;
  Profile_Index *P;
  FILE *kanno, *kdata;

  // Parse arguments
  { int    i, j, k;
    int    flags[128];
    char  *eptr;
    (void) flags;

    ARG_INIT("ShowSeeds");

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
    Trim_DB(db);
    Load_All_Reads(db,0);

    if (verbose)
      fprintf(stderr,"Opening %s.class\n",fk_root);

    // if ((F = gzopen(Catenate(fk_root,"","",".class"), "r")) == NULL)
    //   { fprintf(stderr,"%s: Cannot open %s.class [errno=%d]\n",Prog_Name,fk_root,errno);
    //     exit (1);
    //   }
    // ks = kseq_init(F);

    if (verbose)
      fprintf(stderr,"Opening .%s.class.[anno|data]\n",fk_root);

    kanno = Fopen(Catenate(".",fk_root,"",".class.anno"),"r");
    kdata = Fopen(Catenate(".",fk_root,"",".class.data"),"r");
    if (kanno == NULL || kdata == NULL)
      exit(1);

    int64 cloff;
    if (fseek(kanno,sizeof(int)*2 + sizeof(int64)*db->tfirst,SEEK_SET) < 0)
      { fprintf(stderr,"%s@A: Requies a .class track, run DBclass\n",Prog_Name);
        exit (1);
      }
    if (fread(&cloff,sizeof(int64),1,kanno) != 1)
      { fprintf(stderr,"%s@B: Requies a .class track, run DBclass\n",Prog_Name);
        exit (1);
      }
    fclose(kanno);

    if (fseeko(kdata,cloff,SEEK_SET) < 0)
      { fprintf(stderr,"%s@C: Requies a .class track, run DBclass\n",Prog_Name);
        exit (1);
      }
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
      { rlen = db->reads[id].rlen;

        if (verbose)
          printf("Read %lld (%5d bp):\n",id+1,rlen);

        if (verbose)
          { printf("  POS = ");
            for (int i = 0; i < rlen; i++)
              printf(" %5d",i);
            printf("\n");
          }

        // Prof
        int plen = Fetch_Profile(P,id,db->maxlen+1,profile);
        if (rlen != P->kmer-1+plen)
          { fprintf(stderr,"Length inconsistent\n");
            exit(1);
          }
        if (verbose)
          { printf(" PROF = ");
            for (int i = 0; i < P->kmer-1; i++)
              printf("      ");
            for (int i = 0; i < plen; i++)
              printf(" %5d",profile[i]);
            printf("\n");
          }

        // Class
        // if (kseq_read(ks) < 0)
        //   { fprintf(stderr,"Cannot load seq from %s.class\n",fk_root);
        //     exit(1);
        //   }
        // if (verbose)
        //   { fprintf(stderr,"CLASS = ");
        //     for (int i = 0; i < rlen; i++)
        //       fprintf(stderr," %5c",ks->qual.s[i]);
        //     fprintf(stderr,"\n");
        //   }

        // Seed
        // TODO: Load class without DB? (know clen by anno)
        clen = COMPRESSED_LEN(rlen);
        if (clen <= 0)
          { fprintf(stderr,"clen = 0\n");
            exit(1);
          }
        if (fread(cls,1,clen,kdata) != clen)
          { fprintf(stderr,"%s: Cannot load read %lld\n",Prog_Name,id);
            exit (1);
          }

        Uncompress_Read(rlen,cls);
        c = cls-o;
        if (verbose)
          printf("SEED:\n");
        for (p = o; s[p] != 4; p++)
          if (c[p] > 0)
            { printf("%lld\t%lld\t%c\t%d\t",id+1,p-o,ctos[(int)c[p]],profile[p-o-P->kmer+1]);
              for (int64 pp = p-P->kmer+1; pp <= p; pp++)
                printf("%c",btoc[(int)s[pp]]);
              printf("\n");
            }
        if (p-o != rlen)
          { fprintf(stderr,"rlen (%d) != seq len (%lld)\n",rlen,p-o);
            exit(1);
          }

        o += (rlen+1);

        if (verbose)
          printf("\n");
      }
  }

  fclose(kdata);

  return 0;
}
