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

char *Usage            = "[-vs] [-c<int>] [-r<int(20000)>] [-N<fastk_root>] [-P<tmp_dir(/tmp)>] [-T<int(4)>] "
                         "<source>[.db|.dam|.f[ast][aq][.gz] ...";

int   DEFAULT_NTHREADS = 4;
int   DEFAULT_RLEN     = 20000;
char *DEFAULT_TMP_PATH = "/tmp";

char  stoc[N_STATE]    = {'E','R','H','D'};

char *EXT[N_EXT]       = { ".db", ".dam",
                           ".fastq", ".fasta", ".fq", ".fa",
                           ".fastq.gz", ".fasta.gz", ".fq.gz", ".fa.gz" };

char *osep[N_OTYPE]    = { "/", "/.", "/." };
char *osuf[N_OTYPE]    = { ".class", ".class.data", ".class.anno" };
bool  oann[N_OTYPE]    = { false, false, true };
bool  obin[N_OTYPE]    = { false, true, true };
