#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "ClassPro.h"

static inline cnt_t plus_sigma(cnt_t cnt, int n_sigma)
{ return cnt + (cnt_t)(sqrt(cnt) * n_sigma);
}

static inline cnt_t minus_sigma(cnt_t cnt, int n_sigma)
{
#ifdef DEBUG
  if (cnt <= (cnt_t)(sqrt(cnt) * n_sigma))
    { fprintf(stderr,"cnt (%d) <= its %d sigma\n",cnt,n_sigma);
      exit(1);
    }
#endif
  return cnt - (cnt_t)(sqrt(cnt) * n_sigma);
}

static inline double linear_interpolation(pos_t x, Pos_Cnt pc1, Pos_Cnt pc2)
{
#ifdef DEBUG
  if (!(pc1.pos < x && x < pc2.pos))
    { fprintf(stderr,"Invalid points for interpolation: x1=%d, x=%d, x2=%d\n",pc1.pos,x,pc2.pos);
      exit(1);
    }
#endif
  return (double)pc1.cnt+((double)pc2.cnt-pc1.cnt)*(x-pc1.pos)/(pc2.pos-pc1.pos);
}

static inline double logp_trans(pos_t b, pos_t e, int cb, int ce, cnt_t cov)
{
// #ifdef DEBUG
//   if (!(b < e))
//     { fprintf(stderr,"Violate b (%d) < e (%d)\n",b,e);
//       exit(1);
//     }
// #endif
  return logp_skellam(ce-cb,(double)cov*abs(e-b)/READ_LEN);
}

static inline double p_errorin(enum Etype e, double erate, cnt_t cout, cnt_t cin)
{
#ifdef DEBUG
  if (!(cin <= cout))
    { fprintf(stderr,"Violate cin (%d) <= cout (%d)\n",cin,cout);
      exit(1);
    }
#endif
  return binom_test_g((e == SELF) ? cin : cout-cin,cout,erate,false);
}
