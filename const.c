#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "ClassPro.h"

const char    *Usage            = "[-vs] [-c<int>] [-r<int(20000)>] [-N<fastk_root>] [-P<tmp_dir(./)>] [-T<int(4)>] "
                                  "<source>[.db|.dam|.f[ast][aq][.gz] ...";

const int      DEFAULT_NTHREADS = 4;
const int      DEFAULT_RLEN     = 20000;
const char    *DEFAULT_TMP_PATH = "./";

const int      MERGE_BUF_SIZE   = 4096;   // TODO: larger is faster?

const char     stoc[N_STATE]    = { 'E', 'R', 'H', 'D' };

const char    *EXT[N_EXT]       = { ".db",       ".dam",
                                    ".fastq",    ".fasta",    ".fq",    ".fa",
                                    ".fastq.gz", ".fasta.gz", ".fq.gz", ".fa.gz" };

const Out_Info O_INFO[N_OTYPE]  = { { "/",  ".class",      false, false },
                                    { "/.", ".class.data", false, true  },
                                    { "/.", ".class.anno", true,  true  }  };
