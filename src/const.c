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
                                  "[-P<tmp_dir(./)>] [-N<fastk_root>] [-M<model_path>] "
                                  "<source>[.db|.dam|.f[ast][aq][.gz]";

const char     stoc[N_STATE]    = { 'E', 'R', 'H', 'D' };

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

#define        MAX_KMER_CNT       32767

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

const int      MERGE_BUF_SIZE   = 4096;    // TODO: larger is faster?
const int      MAX_READ_LEN     = 60000;   // For fastx inputs
const int      N_SIGMA_RCOV     = 5;

const int      MAX_N_LC         = 20;      // Max number of bases for a single low-complexity sequence
const int      MAX_N_HC         = 5;       // Max # of bases in a single high-complexity error event
const int      MIN_CNT_CHANGE   = 3;       // Every count change at a wall must be > this
const int      MAX_CNT_CHANGE   = 5;       // Every count change > this becomes a wall candidate(H-cov?)
const double   PE_THRES[N_THRES][N_ETYPE] = { {0.001, 0.05},
                                              {1e-5,  1e-5} };
const double   THRES_DIFF_EO    = -23.025851;   // log(1e-10);
const double   THRES_DIFF_REL   = -9.210340;    // log(1e-4);

const int      OFFSET           = 1000;
const int      N_SIGMA_R        = 2;
const double   R_LOGP           = -10.;
const double   E_PO_BASE        = -10.;
const double   PE_MEAN          = 0.01;
