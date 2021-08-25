#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "ClassPro.h"

#define DIGAMMA_MAX 10000000   // TODO: more efficient way (doubling?)
static double       digamma[DIGAMMA_MAX+1];   // Only for integers

static const int    PMM_MAX_NITER         = 10;
static const int    GAMMA_WEIGHT[2]       = {1, 1};   // Hyperparameters
static const int    ALPHA_WEIGHT          = 1;
static int          a_prior[2], b_prior[2], alpha_prior[2];
int                 lambda_prior[2];   // = [Global H-cov, D-cov]
static double       dg_sum_alpha_prior, eta_weight_k_prior[2], eta_const_k_prior[2];

void precompute_digamma()
{ digamma[1] = -0.57721566490153;
  for (int n = 1; n < DIGAMMA_MAX; n++)
    digamma[n+1] = digamma[n]+(1./(double)n);

  return;
}

void process_global_hist(char *FK_ROOT, int COVERAGE)
{ Histogram *H;
  int64     *hist;

  H = Load_Histogram(FK_ROOT);
  if (H == NULL)
    { fprintf(stderr,"%s: Cannot open %s.hist\n",Prog_Name,FK_ROOT);
      exit(1);
    }
  Modify_Histogram(H,H->low,H->high,0);
  hist = H->hist;

  if (VERBOSE)
    fprintf(stderr,"Global histogram inspection:\n");

  // Find H/D peaks (= mean depths)
  if (COVERAGE != -1)
    { lambda_prior[1] = COVERAGE;
      lambda_prior[0] = COVERAGE >> 1;
      if (VERBOSE)
        fprintf(stderr,"    Specified (H,D) cov   = (%d,%d)\n",lambda_prior[0],lambda_prior[1]);
    }
  else
    { int        maxcnt;               // Tallest peak; Assuming this is H or D
      int64      maxpk;
      int        lmaxcnt, rmaxcnt;     // l for around maxcnt/2, r for around maxcnt*2
      int64      lmaxpk, rmaxpk;       // Assuming taller peak out of l and r is D or H
      int        is_lpeak, is_rpeak;   // Is l/r a peak (rather than monotonously decreasing)?
      double     m, s;

      maxcnt = maxpk = 0;
      for (int i = MAX(2,H->low); i < MIN(1000,H->high); i++)   // NOTE: Assuming coverage at most 1000
        { if (hist[i-1] < hist[i] && hist[i] > hist[i+1] && maxpk < hist[i])
            { maxcnt = i;
              maxpk = hist[i];
            }
        }
      if (maxcnt < 10)
        { fprintf(stderr,"[ERROR] Could not find any peak count >= 10 in the histogram. Revise data and use the `-c` option.");
          exit(1);
        }
      else if (VERBOSE)
        fprintf(stderr,"    Tallest peak count    = %d (# of k-mers = %lld)\n",maxcnt,maxpk);

      m = (double)maxcnt/2;
      s = sqrt(m);
      lmaxcnt = lmaxpk = is_lpeak = 0;
      for (int i = (int)round(m-s); i <= (int)round(m+s); i++)
        { if (lmaxpk < hist[i])
            { lmaxcnt = i;
              lmaxpk = hist[i];
              is_lpeak = (hist[i-1] < hist[i] && hist[i] > hist[i+1]) ? 1 : 0;
            }
        }

      m = (double)maxcnt*2;
      s = sqrt(m);
      rmaxcnt = rmaxpk = is_rpeak = 0;
      for (int i = (int)round(m-s); i <= (int)round(m+s); i++)
        { if (rmaxpk < hist[i])
            { rmaxcnt = i;
              rmaxpk = hist[i];
              is_rpeak = (hist[i-1] < hist[i] && hist[i] > hist[i+1]) ? 1 : 0;
            }
        }

      if (lmaxpk > rmaxpk)   // maxcnt is D-peak and lmaxcnt is H-peak
        { lambda_prior[1] = maxcnt;
          lambda_prior[0] = is_lpeak ? lmaxcnt : (maxcnt >> 1);
        }
      else   // maxcnt is H-peak and rmaxcnt is D-peak
        { lambda_prior[0] = maxcnt;
          lambda_prior[1] = is_rpeak ? rmaxcnt : (maxcnt << 1);
        }

      if (VERBOSE)
        fprintf(stderr,"    Estimated (H,D) cov   = (%d,%d)\n",lambda_prior[0],lambda_prior[1]);
    }

  // Determine hyperparameters
  double totpk[2];   // Estimated total k-mer frequency in H/D components
  double p;

  for (int k = 0; k < 2; k++)
    { b_prior[k] = GAMMA_WEIGHT[k];
      a_prior[k] = lambda_prior[k] * b_prior[k];

      p = exp(logp_poisson(lambda_prior[k],lambda_prior[k]));
      totpk[k] = (double)hist[lambda_prior[k]]/p;
    }
  const int minidx = (totpk[0] < totpk[1]) ? 0 : 1;
  alpha_prior[minidx] = ALPHA_WEIGHT;
  alpha_prior[1-minidx] = (int)(ALPHA_WEIGHT * (totpk[1-minidx] / totpk[minidx]));

  // Precompute values used in VI of PMM
  dg_sum_alpha_prior = digamma[alpha_prior[0] + alpha_prior[1]];
  for (int k = 0; k < 2; k++)
    { eta_weight_k_prior[k] = digamma[a_prior[k]] - log(b_prior[k]);
      eta_const_k_prior[k] = digamma[alpha_prior[k]] - dg_sum_alpha_prior - a_prior[k] / b_prior[k];
    }

  if (VERBOSE)
    { fprintf(stderr,"Hyperparameters etc. for PMM:\n");
      fprintf(stderr,"    lambda_prior          = (%d,%d)\n",lambda_prior[0],lambda_prior[1]);
      fprintf(stderr,"    a_prior               = (%d,%d)\n",a_prior[0],a_prior[1]);
      fprintf(stderr,"    b_prior               = (%d,%d)\n",b_prior[0],b_prior[1]);
      fprintf(stderr,"    alpha_prior           = (%d,%d)\n",alpha_prior[0],alpha_prior[1]);
      fprintf(stderr,"    dg_sum_alpha_prior    = %lf\n",dg_sum_alpha_prior);
      fprintf(stderr,"    eta_weight_k_prior    = (%lf,%lf)\n",eta_weight_k_prior[0],eta_weight_k_prior[1]);
      fprintf(stderr,"    eta_const_k_prior     = (%lf,%lf)\n",eta_const_k_prior[0],eta_const_k_prior[1]);
    }

  Free_Histogram(H);
}
 
int pmm_vi(uint16 *profile, uint16 *nprofile, int plen, double *eta, double lambda[2])
{ int    N;   // # of normal counts, = length of `nprofile` and `eta`/2
  int    is_converged;
  double a[2], b[2], alpha[2], eta_weight_k[2], eta_const_k[2];
  double eta_sum, dg_sum_alpha;

  const int ethres = lambda_prior[0]-3*(int)sqrt(lambda_prior[0]);
  const int rthres = lambda_prior[1]+3*(int)sqrt(lambda_prior[1]);

  for (int k = 0; k < 2; k++)
    { a[k]      = a_prior[k];
      b[k]      = b_prior[k];
      alpha[k]  = alpha_prior[k];
      lambda[k] = lambda_prior[k];
    }

  N = 0;
  for (int i = 0; i < plen; i++)
    if (ethres <= profile[i] && profile[i] <= rthres)
      nprofile[N++] = profile[i];
  if (N < 2)
    return N;

  for (int t = 0; t < PMM_MAX_NITER; t++)
    {
#ifdef DEBUG_PMM
      fprintf(stderr,"[PMM epoch %d]\n",t+1);
      fflush(stderr);
#endif

      // Update eta (assignment weights)
      if (t == 0)
        { dg_sum_alpha = dg_sum_alpha_prior;
          for (int k = 0; k < 2; k++)
            { eta_weight_k[k] = eta_weight_k_prior[k];
              eta_const_k[k] = eta_const_k_prior[k];
            }
        }
      else
        { if (alpha[0]+alpha[1] > DIGAMMA_MAX)
            { fprintf(stderr,"alpha[0]+alpha[1] (%lf) exceeded DIGAMMA_MAX\n",alpha[0]+alpha[1]);
              exit(1);
            }
          dg_sum_alpha = digamma[(int)(alpha[0] + alpha[1])];
          for (int k = 0; k < 2; k++)
            { if (a[k] > DIGAMMA_MAX)
                { fprintf(stderr,"a[%d] (%lf) exceeded DIGAMMA_MAX\n",k,a[k]);
                  exit(1);
                }
              eta_weight_k[k] = digamma[(int)a[k]] - log(b[k]);
              if (alpha[k] > DIGAMMA_MAX)
                { fprintf(stderr,"alpha[%d] (%lf) exceeded DIGAMMA_MAX\n",k,alpha[k]);
                  exit(1);
                }
              eta_const_k[k] = digamma[(int)alpha[k]] - dg_sum_alpha - a[k] / b[k];
            }
        }
      for (int n = 0; n < N; n++)
        { eta_sum = 0.;
          for (int k = 0; k < 2; k++)
            { int idx = (n<<1)|k;
              eta[idx] = exp(nprofile[n] * eta_weight_k[k] + eta_const_k[k]);
              eta_sum += eta[idx];
            }
          for (int k = 0; k < 2; k++)
            eta[(n<<1)|k] /= eta_sum;
        }

/*#ifdef DEBUG_PMM
      for (int n = 0; n < N; n++)
        fprintf(stderr,"n=%d, eta=(%lf,%lf)\n",n,eta[(n<<1)|0],eta[(n<<1)|1]);
#endif*/

      // Update distribution parameters
      // NOTE: Using expected value (as a good proxity for sum of assignments) in alpha updates
      for (int k = 0; k < 2; k++)
        a[k] = b[k] = 0.;
      for (int n = 0; n < N; n++)
        { for (int k = 0; k < 2; k++)
            { a[k] += eta[(n<<1)|k] * nprofile[n];
              b[k] += eta[(n<<1)|k];
            }
        }
      for (int k = 0; k < 2; k++)
        { alpha[k] = b[k] + alpha_prior[k];
          a[k] += a_prior[k];
          b[k] += b_prior[k];
        }

      // Compute lambda (depths)
      is_converged = 1;
      for (int k = 0; k < 2; k++)
        { double _lambda = a[k] / b[k];
          if (fabs(lambda[k] - _lambda) >= 0.1)
            is_converged = 0;
          lambda[k] = _lambda;
        }

#ifdef DEBUG_PMM
      fprintf(stderr,"  alpha  = (%lf,%lf)\n",alpha[0],alpha[1]);
      fprintf(stderr,"  a,b    = [(%lf,%lf), (%lf,%lf)]\n",a[0],b[0],a[1],b[1]);
      fprintf(stderr,"  lambda = (%lf,%lf)\n",lambda[0],lambda[1]);
      fflush(stderr);
#endif

      if (is_converged)
        break;
    }

#if defined(DEBUG_ITER) || defined(DEBUG_PMM)
  if (!is_converged)
    { fprintf(stderr,"PMM Not converged ");
#ifndef DEBUG_ITER
      fprintf(stderr,"\n");
#endif
      fflush(stderr);
    }
#endif

  // Isolate H-cov and D-cov if they are too close
  if (fabs(lambda[0]-lambda[1]) < sqrt(lambda[1]))
    { const double mean = (lambda[0]+lambda[1])/2;
      const double diff_h = fabs(mean-lambda_prior[0]);
      const double diff_d = fabs(mean-lambda_prior[1]);

      if (diff_h < diff_d)
        { lambda[1] = lambda[0]*2;
#if defined(DEBUG_ITER) || defined(DEBUG_PMM)
          fprintf(stderr,"D=H*2 ");
#ifndef DEBUG_ITER
          fprintf(stderr,"\n");
#endif
          fflush(stderr);
#endif
        }
      else
        { lambda[0] = lambda[1]/2;
#if defined(DEBUG_ITER) || defined(DEBUG_PMM)
          fprintf(stderr,"H=D/2 ");
#ifndef DEBUG_ITER
          fprintf(stderr,"\n");
#endif
          fflush(stderr);
#endif
        }
    }

  return N;
}
