#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "ClassPro.h"

/****************************************************************************
 *
 *  ABSOLUTE CONSTANTS
 *
 ****************************************************************************/

const char    *Usage            = "[-vs] [-T<int(4)>] "
                                  "[-c<int>] [-r<int(20000)>] "
                                  "[-P<tmp_dir(./)>] [-N<fastk_root>] "
                                  "<source>[.db|.dam|.f[ast][aq][.gz] ...";

const char     stoc[N_STATE]    = { 'E', 'R', 'H', 'D' };
#ifdef WRITE_TRACK
const char     ctos[128]        = { 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 3, 0, 0, 0,
                                    2, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 1, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0  };
#endif

#define N_EXT 10
const char    *EXT[N_EXT]       = { ".db",       ".dam",
                                    ".fastq",    ".fasta",
                                    ".fq",       ".fa",
                                    ".fastq.gz", ".fasta.gz",
                                    ".fq.gz",    ".fa.gz"     };

const Out_Info O_INFO[N_OTYPE]  = { { "/",  ".class",      false, false },
                                    { "/.", ".class.data", false, true  },
                                    { "/.", ".class.anno", true,  true  }  };

/****************************************************************************
 *
 *  DEFAULT ARGUMENT VALUES
 *
 ****************************************************************************/

const int      DEFAULT_NTHREADS = 4;
const int      DEFAULT_RLEN     = 20000;
const char    *DEFAULT_TMP_PATH = "./";

/****************************************************************************
 *
 *  CONSTANT PARAMETERS
 *
 ****************************************************************************/

const int      MERGE_BUF_SIZE   = 4096;   // TODO: larger is faster?
const int      MAX_READ_LEN     = 50000;   // For fastx inputs
