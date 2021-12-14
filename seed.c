#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "ClassPro.h"

// NOTE: len(profile) == len(class) == plen, len(seq) = plen + K - 1
// NOTE: class[i] in {'E', 'H', 'D', 'R'}
void find_seeds(const char *seq, const uint16 *profile, const char *class, const int plen, const int K)
{ // Calc uniformly sparse count maximizers
  const char C = 'H';
  fprintf(stderr,"\nK=%d\n",K);
  for (int i = 0; i < plen; i++)
    { if (class[i] != C)
        continue;
      fprintf(stderr,"@ %5d: kmer = %.*s, count = %d\n",i,K,seq+i-K+1,profile[i]);
    }

  // Calc uniformly sparse modimizers among the maximizers

  return;
}
