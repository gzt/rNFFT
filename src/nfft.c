
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>  // memset, memcpy
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include<complex.h>

#define NFFT_PRECISION_DOUBLE
#include<fftw3.h>
#include "nfft3mp.h"



#define ALLOC_VECTOR(S, D, ST, DT, C, N)                                       \
  SEXP S;                                                                      \
  PROTECT(S = allocVector(ST, N));                                             \
  DT *D = C(S);

#define ALLOC_MATRIX(S, D, ST, DT, C, NROW, NCOL)                              \
  SEXP S;                                                                      \
  PROTECT(S = allocMatrix(ST, NROW, NCOL));                                    \
  DT *D = C(S);

#define ALLOC_REAL_VECTOR(S, D, N) ALLOC_VECTOR(S, D, REALSXP, double, REAL, N)

#define ALLOC_REAL_MATRIX(S, D, NROW, NCOL)                                    \
  ALLOC_MATRIX(S, D, REALSXP, double, REAL, NROW, NCOL)

#define ALLOC_COMPLEX_VECTOR(S, D, N)                                          \
  ALLOC_VECTOR(S, D, CPLXSXP, Rcomplex, COMPLEX, N)

#define ALLOC_COMPLEX_MATRIX(S, D, NROW, NCOL)                                    \
  ALLOC_MATRIX(S, D, CPLXSXP, Rcomplex, COMPLEX, NROW, NCOL)


void rand_unit_complex(double complex *x, const int n)
{
  int k;

  for (k = 0; k < n; k++)
    x[k] = unif_rand() + I * unif_rand();
}

void rand_shifted_unit_double(double *x, const int n)
{
  int k;

  for (k = 0; k < n; k++)
    x[k] = unif_rand() - (0.5);
}




SEXP test_function(SEXP M, SEXP N){
  NFFT(plan) p;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;
   GetRNGstate();
   NFFT(init_1d)(&p, n, m);
    Rprintf("Initialized. \n");
   rand_shifted_unit_double(p.x, p.M_total);
   Rprintf("Nodes generated. \n");
   for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]);
   nfft_precompute_one_psi(&p); 

   rand_unit_complex(p.f_hat,p.N_total);
  
   for(int i = 0; i < p.N_total; i++) Rprintf("Here the real part of f_hat: %f\n",crealf(p.f_hat[i]));
   
  /** check for valid parameters before calling any trafo/adjoint method */
  error_str = nfft_check(&p);
  if (error_str != 0)
  {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** direct trafo and show the result */
  nfft_trafo_direct(&p);
    
  for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f\n",crealf(p.f[i]));
 
  /** approx. trafo and show the result */
  nfft_trafo(&p);
  
  for(int i = 0; i < p.M_total; i++) Rprintf("Here is nfft %f\n",crealf(p.f[i]));
 
    /** approx. adjoint and show the result */
  nfft_adjoint_direct(&p);
  for(int i = 0; i < p.N_total; i++) Rprintf("Here is adjoint ndft %f\n",crealf(p.f_hat[i]));
 
  /** approx. adjoint and show the result */
  nfft_adjoint(&p);
  for(int i = 0; i < p.N_total; i++) Rprintf("Here is adjoint nfft %f\n",crealf(p.f_hat[i]));
 

  /** finalise the one dimensional plan */
   
    NFFT(finalize)(&p);
   PutRNGstate();
   Rprintf("Finish the program. \n");
   return R_NilValue;
}

