#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "ClassPro.h"

void calc_seq_context(Seq_Ctx *lctx, Seq_Ctx *rctx, char *seq, const int rlen)
{ int       in_hp, in_ds, in_ts;
  const int rlenm1 = rlen-1;

  in_ds = in_ts = 0;
  for (int i = 1; i < rlen; i++)
    { in_hp = (seq[i-1] == seq[i]) ? 1 : 0;
      in_ds = in_ts = 0;

      if (in_hp)
        { lctx[i][HP] = lctx[i-1][HP]+1;
          lctx[i][DS] = rctx[i-1][DS] = 0;
        }
      else
        { lctx[i][HP] = 1;
          lctx[i][DS] = rctx[i-1][DS] = 1;
          for (int j = i-lctx[i-1][HP], n = 0; j < i; j++, n++)
            rctx[j][HP] = lctx[i-1-n][HP];
          if (i >= 3 && seq[i-3] == seq[i-1] && seq[i-2] == seq[i])
            { lctx[i][DS] = lctx[i-2][DS]+1;
              in_ds = 1;
            }
        }

      if (!in_ds)
        { int l = i-1;
          while (lctx[l][DS] > 1)
            l--;
          if (l < i-1)
            for (int j = l-1, n = 0; j < i; j++, n++)
              rctx[j-1][DS] = lctx[i-1-n][DS];
        }

      if (i >= 2)
        { if (in_hp && seq[i-2] == seq[i-1])
            lctx[i][TS] = rctx[i-2][TS] = 0;
          else if (i >= 5 && seq[i-5] == seq[i-2] && seq[i-4] == seq[i-1] && seq[i-3] == seq[i])
            { lctx[i][TS] = lctx[i-3][TS]+1;
              in_ts = 1;
            }
          else
            lctx[i][TS] = rctx[i-1][TS] = rctx[i-2][TS] = 1;

          if (!in_ts)
            { int l = i-1;
              while (lctx[l][TS] > 1)
                l--;
              if (l < i-1)
                for (int j = l-2, n = 0; j < i; j++, n++)
                  rctx[j-2][TS] = lctx[i-1-n][TS];
            }
        }
    }

  for (int j = rlen-lctx[rlenm1][HP], n = 0; j < rlen; j++, n++)
    rctx[j][HP] = lctx[rlenm1-n][HP];

  if (in_ds)
    { int l = rlenm1;
      while (lctx[l][DS] > 1)
        l--;
      if (l < rlenm1)
        for (int j = l-1, n = 0; j < rlen; j++, n++)
          rctx[j-1][DS] = lctx[rlenm1-n][DS];
    }

  if (in_ts)
    { int l = rlenm1;
      while (lctx[l][TS] > 1)
        l--;
      if (l < rlenm1)
        for (int j = l-2, n = 0; j < rlen; j++, n++)
          rctx[j-2][TS] = lctx[rlenm1-n][TS];
    }

  rctx[rlenm1][DS] = rctx[rlenm1][TS] = rctx[rlen-2][TS] = 0;

#ifdef DEBUG_CTX
  const Seq_Ctx *ctx[N_WTYPE]  = {lctx, rctx};
  const char    *dir[N_WTYPE]  = {"L", "R"};
  const char    *name[N_CTYPE] = {"HP", "DS", "TS"};
  const int      W             = 50;

  fprintf(stderr,"Seq  %.*s...%.*s\n",W,seq,W,seq+rlen-W);
  for (int t = HP; t <= TS; t++)
    { for (int d = DROP; d <= GAIN; d++)
        { fprintf(stderr,"%s %s ",name[t],dir[d]);
          for (int i = 0; i < W; i++)
            fprintf(stderr,"%u",ctx[d][i][t]);
          fprintf(stderr,"...");
          for (int i = rlen-W; i < rlen; i++)
            fprintf(stderr,"%u",ctx[d][i][t]);
          fprintf(stderr,"\n");
        }
    }
  fprintf(stderr,"\n");
  fflush(stderr);
#endif

  return;
}
