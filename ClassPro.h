#ifndef _CLASSPRO
#define _CLASSPRO

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// -------- DEBUG FLAGS ---------- //

#define DUP_PROFILE
#undef PARALLEL_WRITE
#define NREAD_PWRITE 100

// #define DEBUG_SINGLE
// #define DEBUG_SINGLE_ID 2
// #undef DEBUG_SMALL
// #define NREAD_SMALL 100

// #define DEBUG
// #define DEBUG_ITER
// #undef DEBUG_BINOM
// #undef DEBUG_CTX
// #undef DEBUG_ERROR
// #undef DEBUG_INTVL
// #undef DEBUG_COR
// #undef DEBUG_PMM
// #undef DEBUG_PROB
// #define DEBUG_REL
// #define DEBUG_UNREL
// #define DEBUG_SLIP
// #undef DEBUG_MERGE

#if defined(DEBUG_SINGLE) || defined(DEBU_SMALL)
#define NO_WRITE
#else
#undef NO_WRITE
#endif

// ------------------------------- //


#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

// Constants
#define MAX_NAME 10000   // .db, .dam header

extern char *Usage;
extern int   DEFAULT_NTHREADS;
extern int   DEFAULT_RLEN;
extern char *DEFAULT_TMP_PATH;
extern int   MERGE_BUF_SIZE;   // for `merge_files`

// Used during multi-thread classification
extern bool  VERBOSE;
extern int   READ_LEN;
extern bool  IS_DB;
extern bool  IS_DAM;

/*******************************************************************************************
 *
 *  CLASSIFICATION
 *
 ********************************************************************************************/

enum State { ERROR, REPEAT, HAPLO, DIPLO, N_STATE };
enum Ctype { HP, DS, TS, N_CTYPE };
enum Etype { SELF, OTHERS, N_ETYPE };
enum Wtype { DROP, GAIN, N_WTYPE };

extern char stoc[N_STATE];

typedef uint8 Seq_Ctx[N_CTYPE];

typedef struct
  { uint8      lmax;    // Maximum feature length considered
    double    *pe;      // Error probability given feature length
    //double    *lpe;     // Log error probability
    //double    *l1mpe;   // Log (1 - error probability)
    double   **pe_bt;
    int    ***cthres;
  } Error_Model;

typedef double P_Error[N_ETYPE][N_WTYPE];

typedef struct
  { int    i;
    int    j;
    double pe;
  } Error_Intvl;

//typedef Error_Intvl Error_Intvls[N_ETYPE];

// TODO: function polymorphism for Intvl and Rel_Intvl?
typedef struct
  { int    i;
    int    j;
    bool   is_rel;   // hard assignment of reliable intervals
    bool   is_err;   // hard assignment of erroneous intervals
    char   asgn;
  } Intvl;

typedef struct
  { int    i;
    int    j;
    int    ci;   // corrected counts
    int    cj;   // @ j-1
    char   asgn;
  } Rel_Intvl;

void calc_seq_context(Seq_Ctx *lctx, Seq_Ctx *rctx, char *seq, const int rlen);

void find_wall(uint16 *profile, int plen, Seq_Ctx *ctx[N_WTYPE],
               Error_Model *emodel, P_Error *perror, P_Error *cerror,
               Error_Intvl *eintvl[N_ETYPE], Intvl *intvl, Rel_Intvl *rintvl,
               int *wall, char *asgn, const int K, int *_N, int *_M);

int pmm_vi(uint16 *profile, uint16 *nprofile, int plen, double *eta, double lambda[2]);

void classify_reliable(Rel_Intvl *rintvl, int M, Intvl *intvl, int N, int plen,
                       P_Error *perror, P_Error *cerror, int hcov, int dcov);

void classify_unreliable(uint16 *profile, int plen, Intvl *intvl, int N,
                         P_Error *perror, int hcov, int dcov);

void remove_slip(uint16 *profile, int plen, Seq_Ctx *ctx[N_WTYPE], char *crack);

/*******************************************************************************************
 *
 *  CLASSIFICATION SETUP
 *
 ********************************************************************************************/

// Arguments for classification (defined for each thread)
typedef struct
  { Profile_Index *P;        // To fetch count profiles
    Error_Model   *emodel;   // Error models for low-complexity bases
    int            beg;      // Reads in [beg,end) are classified in this thread
    int            end;
    gzFile         fxfp;     // FASTX file pointer                                         // TODO: change to void *Seq_Info
    kseq_t        *fxseq;    // To fetch reads from FASTX
    DAZZ_DB       *db;       // To fetch nucleotide sequences from .db/.dam
    DAZZ_STUB     *stub;     // To fetch read names from .db
    FILE          *hdrs;     // To fetch read names from .dam
    FILE          *cfile;    // *.class (fastq-like ascii flie)
    FILE          *dfile;    // .*.class.data (DAZZ_DB track; only when input is .db/.dam)
    FILE          *afile;    // .*.class.anno (DAZZ_DB track; only when input is .db/.dam)
  } Class_Arg;

void precompute_probs();

void process_global_hist(char *FK_ROOT, int COVERAGE);

Error_Model *calc_init_thres();

void free_emodel(Error_Model *emodel);

/*******************************************************************************************
 *
 *  IO
 *
 ********************************************************************************************/

// Supported extension names for input sequence files
#define N_EXT 10
extern char *EXT[N_EXT];

// Output file types and their attributes
enum Otype { CLASS, DATA, ANNO, N_OTYPE };

typedef struct
  { char *sep;
    char *suf;
    bool  is_anno;
    bool  is_bin;
  } Out_Info;

extern Out_Info O_INFO[N_OTYPE];

// Command-line arguments + alpha
typedef struct
  { int    verbose;      // `-v` option
    int    nthreads;     // `-T` option
    int    cov;          // `-c` option
    int    rlen;         // `-r` option
    bool   find_seeds;   // `-s` option
    int    nreads;       // Total number of reads in all input files
    int    nparts;       // Number of reads per thread
    char  *tmp_path;     // `-P` option
    char  *fk_root;      // `-N` option
    char **snames;       // `<source>`; List of input sequence file names
    int    nfiles;       // Length of `snames`
    bool   is_db;        // .db or .dam input?
    bool   is_dam;       // .dam?
  } Arg;

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

#endif // _CLASSPRO
