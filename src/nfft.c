

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <math.h>
#include <string.h>  // memset, memcpy
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include<complex.h>
#include<fftw3.h>

#include "nfft3.h"



void test_function(SEXP M, SEXP N){
   nfft_plan p;
  int m = asInteger(M);
  int n = asInteger(N);
  int i = 0;

   GetRNGstate();
  nfft_init_1d(&p, n, m);
  
 
  for (i = 0; i < p.M_total; i++){
    p.x[i] = unif_rand()-0.5;
    Rprintf("plan.x[i]      : 0x%08x\n", p.x[1]);
  }
  
  /* nfft_precompute_one_psi(&p); */
  

  
  nfft_finalize(&p);
  PutRNGstate();
}

