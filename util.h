#ifndef _UTIL_H
#define _UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "ClassPro.h"

static inline int plus_sigma(int cnt, int n_sigma)
{ return cnt + (int)(n_sigma * sqrt(cnt));
}

static inline int minus_sigma(int cnt, int n_sigma)
{ return cnt - (int)(n_sigma * sqrt(cnt));
}

#endif // _UTIL_H
