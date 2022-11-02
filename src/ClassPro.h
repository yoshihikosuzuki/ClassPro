#ifndef _CLASSPRO_H
#define _CLASSPRO_H

/*******************************************************************************************
 *
 *  DEBUG FLAGS
 *
 ********************************************************************************************/

/* --- Debug modes --- */
// `DEBUG_SINGLE`: Single-read mode. No files are output
#undef DEBUG_SINGLE
const int DEBUG_SINGLE_ID = 11;   // Read ID in single-read mode

/* --- Debug flags --- */
// `DEBUG`: Several assertions used for sanity check
#define DEBUG
#undef DEBUG_ITER
#undef DEBUG_BINOM
#undef DEBUG_EMODEL
#undef DEBUG_CTX
#undef DEBUG_WALL
#undef DEBUG_INTVL
#undef DEBUG_COR
#undef DEBUG_PMM
#undef DEBUG_PROB
#undef DEBUG_REL
#undef DEBUG_UNREL
#undef DEBUG_SLIP

/* --- Developmental or experimental modes --- */
// [DEPRECATED] Experimental mode of parallel write to .class (currently only for Linux & DAZZ_DB inputs)
#undef PARALLEL_WRITE
const int NREAD_PWRITE = 100;   // Number of reads per write in parallel-write mode
// Per-read H/D-cov estimation
#undef DO_PMM

/*******************************************************************************************
 *
 *  HEADERS, TYPES, CONSTANTS
 *
 ********************************************************************************************/

#include "libfastk.h"
#include "DB.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

typedef int    pos_t;
typedef uint16 cnt_t;

enum State { ERROR, REPEAT, HAPLO, DIPLO, N_STATE };
enum Ctype { HP, DS, TS, N_CTYPE };
enum Etype { SELF, OTHERS, N_ETYPE };
enum Wtype { DROP, GAIN, N_WTYPE };

extern const char  stoc[N_STATE];
extern const char  ctos[];
extern const char *Usage;
extern const int   DEFAULT_NTHREADS;
extern const int   DEFAULT_RLEN;
extern const char *DEFAULT_TMP_PATH;
extern const int   MERGE_BUF_SIZE;
extern const int   MAX_READ_LEN;
extern const int   N_SIGMA_RCOV;

/*******************************************************************************************
 *
 *  SHARED PARAMETERS
 *
 ********************************************************************************************/

extern bool  VERBOSE;
extern int   READ_LEN;
extern bool  IS_DB;
extern bool  IS_DAM;
extern cnt_t GLOBAL_COV[N_STATE];

/*******************************************************************************************
 *
 *  K-mer count histogram manipulations (hist.c)
 *
 ********************************************************************************************/

extern int lambda_prior[2];   // Global H-cov & D-cov estimated from global count histogram

void precompute_digamma();

void process_global_hist(char *FK_ROOT, int COVERAGE);

typedef struct
  { cnt_t  *nprofile;   // Normal counts for PMM
    double *eta;
  } PMM_Arg;

PMM_Arg *alloc_pmm_arg(int rlen_max);
void free_pmm_arg(PMM_Arg *arg);

int pmm_vi(PMM_Arg *arg, cnt_t *profile, int plen, double lambda[2]);

/*******************************************************************************************
 *
 *  Sequence context (context.c)
 *
 ********************************************************************************************/

typedef uint8 Seq_Ctx[N_CTYPE];

void calc_seq_context(Seq_Ctx *lctx, Seq_Ctx *rctx, char *seq, const int rlen);

/*******************************************************************************************
 *
 *  Wall detection with adaptive error model (wall.c)
 *
 ********************************************************************************************/

enum ThresT { INIT, FINAL, N_THRES };

extern const int    MAX_N_LC;
extern const double PE_THRES[N_THRES][N_ETYPE];

// extern uint8 CMAX;   // Max cout for precomputation of count thresholds

typedef struct
  { uint8      lmax;   // Maximum feature length considered
    double    *pe;         // Error probability given feature length
    uint8  ****cthres;     // Threshold of count change; emodel[ctype].cthres[l][cout][thresT][etype] = cin threshold  // TODO: better order? define outside Error_Model?
  } Error_Model;

Error_Model *calc_init_thres(const char *name);
void free_emodel(Error_Model *emodel);

extern const int    MAX_N_HC;
extern const int    MIN_CNT_CHANGE;
extern const int    MAX_CNT_CHANGE;
extern const double THRES_DIFF_EO;
extern const double THRES_DIFF_REL;

typedef double P_Error[N_ETYPE][N_WTYPE];

//typedef Error_Intvl Error_Intvls[N_ETYPE];

typedef struct
  { double b;
    double e;
  } Perror_O;

typedef struct
  { pos_t  b;
    pos_t  e;
    double pe;
  } Error_Intvl;

typedef struct
  { pos_t    b;
    pos_t    e;
    cnt_t    cb;
    cnt_t    ce;
    cnt_t    ccb;
    cnt_t    cce;
    bool     is_rel;   // hard assignment of reliable intervals
    double   pe;     // TODO: logp?
    Perror_O pe_o;
    char     asgn;
  } Intvl;

typedef struct
  { char        *wall;
    Error_Intvl *eintvl;
    Error_Intvl *ointvl;
    P_Error     *perror;   // perror[i][etype][wtype]
  } Wall_Arg;

Wall_Arg *alloc_wall_arg(int rlen_max);
void free_wall_arg(Wall_Arg *arg);

int find_wall(Wall_Arg *arg, Intvl *intvl, cnt_t *profile, int plen,
               Seq_Ctx *ctx[N_WTYPE], Error_Model *emodel, int K);

int find_rel_intvl(Intvl *intvl, int N, Intvl *rintvl, cnt_t *profile, Seq_Ctx *ctx[N_WTYPE], int K);

/*******************************************************************************************
 *
 *  Interval classification (class_rel.c, class_unrel.c)
 *
 ********************************************************************************************/

extern const int    OFFSET;
extern const int    N_SIGMA_R;
extern const double R_LOGP;
extern const double E_PO_BASE;
extern const double PE_MEAN;

extern double DR_RATIO;

typedef struct
  { pos_t pos;
    cnt_t cnt;
  } Pos_Cnt;

typedef Pos_Cnt CovsT[N_STATE];

#define REL_IDX(i, s) (i)*N_STATE+(s)

typedef struct
  { bool     FORWARD;
    cnt_t   *COV;
    double  *dp;
    CovsT   *st;
    char   **bt;
    double  *dh_ratio;
    bool    *rpos;
    Intvl   *intvl;
  } Rel_Arg;

Rel_Arg *alloc_rel_arg(int rlen_max);
void free_rel_arg(Rel_Arg *arg, int rlen_max);

void classify_rel(Rel_Arg *arg, Intvl *rintvl, int M, Intvl *intvl, int N, int plen);
void classify_unrel(Intvl *intvl, int N);
// void remove_slip(uint16 *profile, int plen, Seq_Ctx *ctx[N_WTYPE], char *crack);

/*******************************************************************************************
 *
 *  Seed selection (seed.c)
 *
 ********************************************************************************************/

// Compressed profile segment
typedef struct {
  int b;     // [b..e) has tie counts of a specific class, or does not belong to the class
  int e;
  int cnt;   // Tie count of this segment if valid, otherwise -1
  int nw;    // Number of segment windows in which this segment is MM
  bool is_seed;
} seg_t;

// Segment for deque
typedef struct {
  int seg_id;   // index of a segment
  int b;        // begin position of the segment
  int e;        // end position, used for `last_oor_pos`
  int cnt;      // (tie) count
} hmer_t;

#include "kdq.h"
KDQ_INIT(hmer_t);

// Simple interval for masked region on a profile
typedef struct {
  int b;
  int e;
} intvl_t;

void find_seeds(kdq_t(hmer_t) *Q,
                const char    *_seq,
                const char    *class,
                const uint16  *profile,
                seg_t         *cprofile,
                int           *hash,
                int           *sasgn,
                intvl_t       *mintvl,
                const int      plen,
                const int      K,
                FILE          *ranno,
                FILE          *rdata,
                int64         *ridx);

/*******************************************************************************************
 *
 *  Main routine for classification (ClassPro.c)
 *
 ********************************************************************************************/

// Command-line arguments + alpha
typedef struct
  { int    verbose;      // `-v` option
    int    nthreads;     // `-T` option
    int    rlen;         // `-r` option
    cnt_t  cov;          // `-c` option
    int    nreads;       // Total number of reads in all input files
    int    nparts;       // Number of reads per thread
    int    nfiles;       // Length of `snames`
    char **snames;       // `<source>`; List of input sequence file names
    char  *tmp_path;     // `-P` option
    char  *fk_root;      // `-N` option
    char  *out_root;     // `out_root`.class is the output file
    char  *model_path;   // `-M` option
    bool   is_db;        // .db or .dam input?
    bool   is_dam;       // .dam?
    bool   find_seeds;   // `-s` option
  } Arg;

// Arguments for classification (defined for each thread)
typedef struct
  { Profile_Index *P;        // To fetch count profiles
    Error_Model   *emodel;   // Error models for low-complexity bases
    int            beg;      // Reads in [beg,end) are classified in this thread
    int            end;
    // Below are for IO
    gzFile         fxfp;     // FASTX file pointer
    kseq_t        *fxseq;    // To fetch reads from FASTX
    DAZZ_DB       *db;       // To fetch nucleotide sequences from .db/.dam
    DAZZ_STUB     *stub;     // To fetch read names from .db
    FILE          *hdrs;     // To fetch read names from .dam
    FILE          *cfile;    // *.class (fastq-like ascii flie)
    FILE          *dfile;    // .*.class.data (DAZZ_DB track; only when input is .db/.dam)
    FILE          *afile;    // .*.class.anno (DAZZ_DB track; only when input is .db/.dam)
    FILE          *rdata;    // .*.rep.data (DAZZ_DB mask track)
    FILE          *ranno;    // .*.rep.anno
  } Class_Arg;

/*******************************************************************************************
 *
 *  IO (io.c)
 *
 ********************************************************************************************/

// Supported extension names for input sequence files
#define        N_EXT              10
const char    *EXT[N_EXT]       = { ".db",       ".dam",
                                    ".fastq",    ".fasta",
                                    ".fq",       ".fa",
                                    ".fastq.gz", ".fasta.gz",
                                    ".fq.gz",    ".fa.gz"     };

// Output file types and their attributes
enum Otype { CLASS, DATA, ANNO, RDATA, RANNO, N_OTYPE };

typedef struct
  { char *sep;
    char *suf;
    bool  is_anno;   // DAZZ_DB data track?
    // bool  is_mask;   // DAZZ_DB mask track?
    bool  is_bin;    // binary file?
  } Out_Info;

const Out_Info O_INFO[N_OTYPE]  = { { "/",  ".class",      false, false },
                                    { "/.", ".class.data", false, true  },
                                    { "/.", ".class.anno", true,  true  },
                                    { "/.", ".rep.data",   false, true  },
                                    { "/.", ".rep.anno",   true,  true  }  };

// const Out_Info O_INFO[N_OTYPE]  = { { "/",  ".class",      false, false, false },
//                                     { "/.", ".class.data", false, false, true  },
//                                     { "/.", ".class.anno", true,  false, true  },
//                                     { "/.", ".rep.data",   false, true,  true  },
//                                     { "/.", ".rep.anno",   true,  true,  true  }  };

// Arguments for merging intermediate output files (defined for each `Otype`)
typedef struct
  { char **onames;   // List of intermediate output file names per thread
    char  *ofinal;   // Name of final output file
    int    nfiles;   // Number of intermediate output files
    bool   is_bin;   // Binary file?
  } Merge_Arg;

void prepare_param(Arg *arg, Class_Arg *paramc, Merge_Arg *paramm);
void free_param(Arg *arg, Class_Arg *paramc, Merge_Arg *paramm);

void *merge_anno(void *arg);
void *merge_files(void *arg);

#endif // _CLASSPRO_H
