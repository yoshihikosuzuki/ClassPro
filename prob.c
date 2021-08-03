#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "prob.h"
#include "bessel.c"

double logfact[32768];

void precompute_logfact()
{ for (int n = 1; n < 32768; n++)
    logfact[n] = logfact[n-1]+log(n);

  return;
}
