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

char *Usage         = "[-vs] [-T<int(4)>] [-c<int>] [-r<int(20000)>] [-N<fastk_root>] "
                      "<source>[.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz] ...";

int DEFAULT_RLEN    = 20000;

char stoc[N_STATE]  = {'E','R','H','D'};

char *EXT[N_EXT]    = { ".db",
                        ".fastq", ".fasta", ".fq", ".fa",
                        ".fastq.gz", ".fasta.gz", ".fq.gz", ".fa.gz" };

char *osep[N_OTYPE] = { "/", "/.", "/." };
char *osuf[N_OTYPE] = { ".class", ".class.data", ".class.anno" };
bool  obin[N_OTYPE] = { false, true, true };
