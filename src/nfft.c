
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

   GetRNGstate();
   NFFT(init_1d)(&p, n, m);
    Rprintf("Got here\n");
   rand_shifted_unit_double(p.x, p.M_total);
   Rprintf("didnt' get here\n");
   for(int i = 0; i < p.M_total; i++) Rprintf("Here is %f\n",p.x[i]);
   nfft_precompute_one_psi(&p); 
   

   
   //  NFFT(finalize)(&p);
   PutRNGstate();
   Rprintf("Do I get here?\n");
   return R_NilValue;
}

