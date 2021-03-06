
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memset, memcpy

#define NFFT_PRECISION_DOUBLE
#include <fftw3.h>
#include <nfft3.h>
#include <nfft3mp.h>

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

#define ALLOC_COMPLEX_MATRIX(S, D, NROW, NCOL)                                 \
  ALLOC_MATRIX(S, D, CPLXSXP, Rcomplex, COMPLEX, NROW, NCOL)

void rand_unit_complex(double complex *x, const int n) {
  int k;

  for (k = 0; k < n; k++)
    x[k] = unif_rand() + I * unif_rand();
}

void rand_shifted_unit_double(double *x, const int n) {
  int k;

  for (k = 0; k < n; k++)
    x[k] = unif_rand() - (0.5);
}

SEXP rndft_1d(SEXP X, SEXP FHAT, SEXP M, SEXP N) {

  NFFT(plan) p;
  int j;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;
  NFFT(init_1d)(&p, n, m);

  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < m; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < m; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    error("'X' must be real or complex.");
  }
  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfft_precompute_one_psi(&p);

  Rcomplex *ffhat = COMPLEX(FHAT);
  for (j = 0; j < n; j++) {
    p.f_hat[j] = ffhat[j].r + I * ffhat[j].i;
  }
  // for(int i = 0; i < p.N_total; i++) Rprintf("Here is f_hat: %f + %f
  // i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFFT(trafo_direct)(&p);

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f + %f
  // i\n",crealf(p.f[i]),cimagf(p.f[i]));

  ALLOC_COMPLEX_VECTOR(F, ret, m);
  for (i = 0; i < m; ++i) {
    ret[i].r = creal(p.f[i]);
    ret[i].i = cimag(p.f[i]);
  }
  NFFT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfft_1d(SEXP X, SEXP FHAT, SEXP M, SEXP N) {

  NFFT(plan) p;
  int j;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;
  NFFT(init_1d)(&p, n, m);

  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < m; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.M_total; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    error("'X' must be real or complex.");
  }
  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfft_precompute_one_psi(&p);

  Rcomplex *ffhat = COMPLEX(FHAT);
  for (j = 0; j < p.N_total; j++) {
    p.f_hat[j] = ffhat[j].r + I * ffhat[j].i;
  }
  // for(int i = 0; i < p.N_total; i++) Rprintf("Here is f_hat: %f + %f
  // i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFFT(trafo)(&p);

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f + %f
  // i\n",crealf(p.f[i]),cimagf(p.f[i]));

  ALLOC_COMPLEX_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i].r = creal(p.f[i]);
    ret[i].i = cimag(p.f[i]);
  }
  NFFT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfft_adjoint_1d(SEXP X, SEXP F, SEXP M, SEXP N) {

  NFFT(plan) p;
  int j;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;
  NFFT(init_1d)(&p, n, m);

  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < m; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < m; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    error("'X' must be real or complex.");
  }
  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfft_precompute_one_psi(&p);

  Rcomplex *ff = COMPLEX(F);
  for (j = 0; j < m; j++) {
    p.f[j] = ff[j].r + I * ff[j].i;
  }
  // for(int i = 0; i < p.N_total; i++) Rprintf("Here is f_hat: %f + %f
  // i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /**  adjoint trafo and show the result */
  NFFT(adjoint)(&p);

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f + %f
  // i\n",crealf(p.f[i]),cimagf(p.f[i]));

  ALLOC_COMPLEX_VECTOR(FHAT, ret, n);
  for (i = 0; i < n; ++i) {
    ret[i].r = creal(p.f_hat[i]);
    ret[i].i = cimag(p.f_hat[i]);
  }
  NFFT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return FHAT;
}

SEXP rndft_adjoint_1d(SEXP X, SEXP F, SEXP M, SEXP N) {

  NFFT(plan) p;
  int j;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;
  NFFT(init_1d)(&p, n, m);

  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < m; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < m; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    error("'X' must be real or complex.");
  }
  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfft_precompute_one_psi(&p);

  Rcomplex *ff = COMPLEX(F);
  for (j = 0; j < m; j++) {
    p.f[j] = ff[j].r + I * ff[j].i;
  }
  // for(int i = 0; i < p.N_total; i++) Rprintf("Here is f_hat: %f + %f
  // i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /**  adjoint trafo and show the result */
  NFFT(adjoint_direct)(&p);

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f + %f
  // i\n",crealf(p.f[i]),cimagf(p.f[i]));

  ALLOC_COMPLEX_VECTOR(FHAT, ret, n);
  for (i = 0; i < n; ++i) {
    ret[i].r = creal(p.f_hat[i]);
    ret[i].i = cimag(p.f_hat[i]);
  }
  NFFT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return FHAT;
}

SEXP rnfft_solver_1d(SEXP X, SEXP Y, SEXP M, SEXP N, SEXP eps,
                     SEXP iterations) {
  /* Y must be length M  */
  NFFT(plan) p;
  SOLVER(plan_complex) ip; /**< plan for the inverse nfft       */
  int k, l;
  int m = asInteger(M);
  int n = asInteger(N);
  int iter = asInteger(iterations);
  double epsilon = REAL(eps)[0];
  const char *error_str;

  NFFT(init_1d)(&p, n, m);

  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < m; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (i = 0; i < m; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    Rf_error("'X' must be real or complex.");
  }

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfft_precompute_one_psi(&p);
  /** initialise inverse plan */
  SOLVER(init_complex)(&ip, (NFFT(mv_plan_complex) *)(&p));

  if (CPLXSXP == TYPEOF(Y)) {
    Rcomplex *yy = COMPLEX(Y);
    for (i = 0; i < p.M_total; ++i) {
      ip.y[i] = yy[i].r + I * yy[i].i;
    }
  } else if (REALSXP == TYPEOF(Y)) {
    double *yy = REAL(X);
    for (int i = 0; i < p.M_total; ++i) {
      ip.y[i] = yy[i] + I * 0;
    }
  } else {
    error("'Y' must be real or complex.");
  }
  /** init pseudo random samples and show them */
  /* rand_unit_complex(ip.y, p.M_total); */
  /* NFFT(vpr_complex)(ip.y, p.M_total, "Given data, vector y"); */

  /** initialise some guess f_hat_0 and solve */
  for (k = 0; k < p.N_total; k++)
    ip.f_hat_iter[k] = NFFT_K(0.0);

  /* NFFT(vpr_complex)(ip.f_hat_iter, p.N_total,      "Initial guess, vector
   * f_hat_iter"); */

  /* for(int i = 0; i < p.N_total; i++) Rprintf("Initial guess: %f + %f i\n", */
  /* 					     crealf(ip.f_hat_iter[i]),
   * cimagf(ip.f_hat_iter[i])); */

  /** check for valid parameters before calling any trafo/adjoint method */
  error_str = NFFT(check)(&p);
  if (error_str != 0) {
    printf("Error in nfft module: %s\n", error_str);
    return R_NilValue;
  }

  NFFT_CSWAP(ip.f_hat_iter, p.f_hat);
  NFFT(trafo)(&p);
  /* NFFT(vpr_complex)(p.f, p.M_total, "Data fit, vector f"); */
  /* for(int i = 0; i < p.M_total; i++) Rprintf("Initial guess: %f + %f i\n", */
  /* 					     crealf(p.f[i]), cimagf(p.f[i])); */

  NFFT_CSWAP(ip.f_hat_iter, p.f_hat);

  SOLVER(before_loop_complex)(&ip);
  /* printf("\n Residual r=%" NFFT__FES__ "\n", ip.dot_r_iter); */

  for (l = 0; ip.dot_r_iter > epsilon && l < iter; l++) {
    /* printf("\n********** Iteration l=%d **********\n", l); */
    SOLVER(loop_one_step_complex)(&ip);
    /* NFFT(vpr_complex)(ip.f_hat_iter, p.N_total, */
    /*     "Approximate solution, vector f_hat_iter"); */
    /* for(int i = 0; i < p.N_total; i++) Rprintf("Approx solution: %f + %f
     * i\n", */
    /* crealf(ip.f_hat_iter[i]), cimagf(ip.f_hat_iter[i])); */
    NFFT_CSWAP(ip.f_hat_iter, p.f_hat);
    NFFT(trafo)(&p);
    /* NFFT(vpr_complex)(p.f, p.M_total, "Data fit, vector f"); */
    /* for(int i = 0; i < p.M_total; i++) Rprintf("Data fit, vector f: %f + %f
     * i\n", */
    /* crealf(p.f[i]), cimagf(p.f[i])); */
    NFFT_CSWAP(ip.f_hat_iter, p.f_hat);

    /* printf("\n Residual r=%"  NFFT__FES__ "\n", ip.dot_r_iter); */
  }

  ALLOC_COMPLEX_VECTOR(FHAT, ret, p.N_total);
  for (i = 0; i < n; ++i) {
    ret[i].r = creal(ip.f_hat_iter[i]);
    ret[i].i = cimag(ip.f_hat_iter[i]);
  }

  SOLVER(finalize_complex)(&ip);
  NFFT(finalize)(&p);
  UNPROTECT(1); /* s_ret */
  return FHAT;
}

SEXP rnfct_1d(SEXP X, SEXP FHAT, SEXP M, SEXP N) {

  NFCT(plan) p;
  int j;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;
  NFCT(init_1d)(&p, n, m);

  int i = 0;
  if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.M_total; ++i) {
      p.x[i] = xx[i];
    }
  } else {
    error("'X' must be real.");
  }
  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfct_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < n; j++) {
    p.f_hat[j] = ffhat[j];
  }
  // for(int i = 0; i < p.N_total; i++) Rprintf("Here is f_hat: %f + %f
  // i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  error_str = nfct_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfct module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFCT(trafo)(&p);

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f + %f
  // i\n",crealf(p.f[i]),cimagf(p.f[i]));

  ALLOC_REAL_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i] = (p.f[i]);
  }
  NFCT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfct_adjoint_1d(SEXP X, SEXP F, SEXP M, SEXP N) {

  NFCT(plan) p;
  int j;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;

  NFCT(init_1d)(&p, n, m);

  int i = 0;
  if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.M_total; ++i) {
      p.x[i] = xx[i];
    }
  } else {
    error("'X' must be real.");
  }
  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfct_precompute_one_psi(&p);

  double *ff = REAL(F);
  for (j = 0; j < p.M_total; j++) {
    p.f[j] = ff[j];
  }
  // for(int i = 0; i < p.N_total; i++) Rprintf("Here is f_hat: %f + %f
  // i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  error_str = nfct_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfct module,\n");
    return R_NilValue;
  }

  /**  adjoint trafo and show the result */
  NFCT(adjoint)(&p);

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f + %f
  // i\n",crealf(p.f[i]),cimagf(p.f[i]));

  ALLOC_REAL_VECTOR(FHAT, ret, p.N_total);
  for (i = 0; i < p.N_total; ++i) {
    ret[i] = (p.f_hat[i]);
  }
  NFCT(finalize)(&p);
  UNPROTECT(1); /* s_ret */
  return FHAT;
}

SEXP rnfst_1d(SEXP X, SEXP FHAT, SEXP M, SEXP N) {

  NFST(plan) p;
  int j;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;

  NFST(init_1d)(&p, n, m);

  int i = 0;
  if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.M_total; ++i) {
      p.x[i] = xx[i];
    }
  } else {
    error("'X' must be real.");
  }
  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfst_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.N_total; j++) {
    p.f_hat[j] = ffhat[j];
  }
  // for(int i = 0; i < p.N_total; i++) Rprintf("Here is f_hat: %f + %f
  // i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  error_str = nfst_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfst module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFST(trafo)(&p);

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f + %f
  // i\n",crealf(p.f[i]),cimagf(p.f[i]));

  ALLOC_REAL_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i] = (p.f[i]);
  }
  NFST(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfst_adjoint_1d(SEXP X, SEXP F, SEXP M, SEXP N) {

  NFST(plan) p;
  int j;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;
  NFST(init_1d)(&p, n, m);

  int i = 0;
  if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.M_total; ++i) {
      p.x[i] = xx[i];
    }
  } else {
    error("'X' must be real.");
  }
  // for(int i = 0; i < p.M_total; i++) Rprintf("Here are the x[i]:
  // %f\n",p.x[i]);
  nfst_precompute_one_psi(&p);

  double *ff = REAL(F);
  for (j = 0; j < p.M_total; j++) {
    p.f[j] = ff[j];
  }
  // for(int i = 0; i < p.N_total; i++) Rprintf("Here is f_hat: %f + %f
  // i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  error_str = nfst_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfst module,\n");
    return R_NilValue;
  }

  /**  adjoint trafo and show the result */
  NFST(adjoint)(&p);

  // for(int i = 0; i < p.M_total; i++) Rprintf("Here is ndft %f + %f
  // i\n",crealf(p.f[i]),cimagf(p.f[i]));

  ALLOC_REAL_VECTOR(FHAT, ret, p.N_total);
  for (i = 0; i < p.N_total; ++i) {
    ret[i] = (p.f_hat[i]);
  }
  NFST(finalize)(&p);
  UNPROTECT(1); /* s_ret */
  return FHAT;
}

SEXP test_function(SEXP M, SEXP N) {
  NFFT(plan) p;
  int m = asInteger(M);
  int n = asInteger(N);
  const char *error_str;
  GetRNGstate();
  NFFT(init_1d)(&p, n, m);
  Rprintf("Initialized. \n");
  rand_shifted_unit_double(p.x, p.M_total);
  Rprintf("Nodes generated. \n");
  for (int i = 0; i < p.M_total; i++)
    Rprintf("Here are the x[i]: %f\n", p.x[i]);
  nfft_precompute_one_psi(&p);

  rand_unit_complex(p.f_hat, p.N_total);

  for (int i = 0; i < p.N_total; i++)
    Rprintf("Here is f_hat: %f + %f i\n", crealf(p.f_hat[i]),
            cimagf(p.f_hat[i]));

  /** check for valid parameters before calling any trafo/adjoint method */
  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** direct trafo and show the result */
  nfft_trafo_direct(&p);

  for (int i = 0; i < p.M_total; i++)
    Rprintf("Here is ndft %f + %f i\n", crealf(p.f[i]), cimagf(p.f[i]));

  /** approx. trafo and show the result */
  nfft_trafo(&p);

  for (int i = 0; i < p.M_total; i++)
    Rprintf("Here is nfft  %f + %f i\n", crealf(p.f[i]), cimagf(p.f[i]));

  /** approx. adjoint and show the result */
  nfft_adjoint_direct(&p);
  for (int i = 0; i < p.N_total; i++)
    Rprintf("Here is adjoint ndft %f + %f i\n", crealf(p.f_hat[i]),
            cimagf(p.f_hat[i]));

  /** approx. adjoint and show the result */
  nfft_adjoint(&p);
  for (int i = 0; i < p.N_total; i++)
    Rprintf("Here is adjoint nfft %f + %f i\n", crealf(p.f_hat[i]),
            cimagf(p.f_hat[i]));

  /** finalise the one dimensional plan */

  NFFT(finalize)(&p);
  PutRNGstate();
  Rprintf("Finish the program. \n");

  return R_NilValue;
}

SEXP solvetest(SEXP m, SEXP n, SEXP iterations) {
  int iter = asInteger(iterations);
  int M = asInteger(m);
  int N = asInteger(n);

  int k, l;                /**< index for nodes, freqencies,iter*/
  NFFT(plan) p;            /**< plan for the nfft               */
  SOLVER(plan_complex) ip; /**< plan for the inverse nfft       */
  const char *error_str;

  /** initialise an one dimensional plan */
  NFFT(init_1d)(&p, N, M);
  GetRNGstate();
  /** init pseudo random nodes */
  rand_shifted_unit_double(p.x, p.M_total);

  /** precompute psi, the entries of the matrix B */

  NFFT(precompute_one_psi)(&p);

  /** initialise inverse plan */
  SOLVER(init_complex)(&ip, (NFFT(mv_plan_complex) *)(&p));

  /** init pseudo random samples and show them */
  rand_unit_complex(ip.y, p.M_total);
  /* NFFT(vpr_complex)(ip.y, p.M_total, "Given data, vector y"); */
  PutRNGstate();
  for (int i = 0; i < p.M_total; i++)
    Rprintf("Here is data vector y: %f + %f i\n", crealf(ip.y[i]),
            cimagf(ip.y[i]));

  /** initialise some guess f_hat_0 and solve */
  for (k = 0; k < p.N_total; k++)
    ip.f_hat_iter[k] = NFFT_K(0.0);

  /* NFFT(vpr_complex)(ip.f_hat_iter, p.N_total,      "Initial guess, vector
   * f_hat_iter"); */

  for (int i = 0; i < p.N_total; i++)
    Rprintf("Initial guess: %f + %f i\n", crealf(ip.f_hat_iter[i]),
            cimagf(ip.f_hat_iter[i]));

  /** check for valid parameters before calling any trafo/adjoint method */
  error_str = NFFT(check)(&p);
  if (error_str != 0) {
    printf("Error in nfft module: %s\n", error_str);
    return R_NilValue;
  }

  NFFT_CSWAP(ip.f_hat_iter, p.f_hat);
  NFFT(trafo)(&p);
  /* NFFT(vpr_complex)(p.f, p.M_total, "Data fit, vector f"); */
  for (int i = 0; i < p.M_total; i++)
    Rprintf("Initial guess: %f + %f i\n", crealf(p.f[i]), cimagf(p.f[i]));

  NFFT_CSWAP(ip.f_hat_iter, p.f_hat);

  SOLVER(before_loop_complex)(&ip);
  printf("\n Residual r=%" NFFT__FES__ "\n", ip.dot_r_iter);

  for (l = 0; l < iter; l++) {
    printf("\n********** Iteration l=%d **********\n", l);
    SOLVER(loop_one_step_complex)(&ip);
    /* NFFT(vpr_complex)(ip.f_hat_iter, p.N_total, */
    /*     "Approximate solution, vector f_hat_iter"); */
    for (int i = 0; i < p.N_total; i++)
      Rprintf("Approx solution: %f + %f i\n", crealf(ip.f_hat_iter[i]),
              cimagf(ip.f_hat_iter[i]));
    NFFT_CSWAP(ip.f_hat_iter, p.f_hat);
    NFFT(trafo)(&p);
    /* NFFT(vpr_complex)(p.f, p.M_total, "Data fit, vector f"); */
    for (int i = 0; i < p.M_total; i++)
      Rprintf("Data fit, vector f: %f + %f i\n", crealf(p.f[i]),
              cimagf(p.f[i]));
    NFFT_CSWAP(ip.f_hat_iter, p.f_hat);

    printf("\n Residual r=%" NFFT__FES__ "\n", ip.dot_r_iter);
  }

  SOLVER(finalize_complex)(&ip);
  NFFT(finalize)(&p);

  Rprintf("Finish the program. \n");
  return R_NilValue;
}

void nfft_2dtest() {
  int N[2], n[2], M;

  NFFT(plan) p;

  const char *error_str;

  N[0] = 32;
  n[0] = 64;
  N[1] = 14;
  n[1] = 32;
  M = N[0] * N[1];

  GetRNGstate();

  /** init a two dimensional plan */
  NFFT(init_guru)
  (&p, 2, N, M, n, 7,
   PRE_PHI_HUT | PRE_FULL_PSI | MALLOC_F_HAT | MALLOC_X | MALLOC_F | FFTW_INIT |
       FFT_OUT_OF_PLACE,
   FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  /* NFFT(init_2d)(&p, N[0], N[1], M); */

  /** init pseudo random nodes */
  rand_shifted_unit_double(p.x, p.d * p.M_total);
  /* NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total); */
  Rprintf("Here is M_total and N_total: %d, %d\n", p.M_total, p.N_total);
  /** precompute psi, the entries of the matrix B */
  /* if(p.flags & PRE_ONE_PSI) */
  NFFT(precompute_one_psi)(&p);

  /** init pseudo random Fourier coefficients and show them */
  /* NFFT(vrand_unit_complex)(p.f_hat, p.N_total); */
  rand_unit_complex(p.f_hat, p.N_total);
  for (int i = 0; i < 7; i++)
    Rprintf("Here is x: %f\n", (p.x[i]));
  for (int i = 0; i < 7; i++)
    Rprintf("Here is f_hat: %f + %f i\n", crealf(p.f_hat[i]),
            cimagf(p.f_hat[i]));
  /* NFFT(vpr_complex)(p.f_hat,K, "given Fourier coefficients, vector f_hat
   * (first few entries)"); */

  /** check for valid parameters before calling any trafo/adjoint method */
  error_str = NFFT(check)(&p);
  if (error_str != 0) {
    printf("Error in nfft module: %s\n", error_str);
    return;
  }

  /** direct trafo and show the result */

  NFFT(trafo_direct)(&p);

  for (int i = 0; i < 7; i++)
    Rprintf("Here is ndft vector f: %f + %f i\n", crealf(p.f[i]),
            cimagf(p.f[i]));

  /* NFFT(vpr_complex)(p.f, K, "ndft, vector f (first few entries)"); */

  /** approx. trafo and show the result */

  NFFT(trafo)(&p);

  /* NFFT(vpr_complex)(p.f, K, "nfft, vector f (first few entries)"); */
  for (int i = 0; i < 7; i++)
    Rprintf("Here is nfft vector f: %f + %f i\n", crealf(p.f[i]),
            cimagf(p.f[i]));

  /** direct adjoint and show the result */

  NFFT(adjoint_direct)(&p);
  for (int i = 0; i < 7; i++)
    Rprintf("Here is adjoint ndft vector f_hat: %f + %f i\n",
            crealf(p.f_hat[i]), cimagf(p.f_hat[i]));

  /* NFFT(vpr_complex)(p.f_hat, K, "adjoint ndft, vector f_hat (first few
   * entries)"); */

  /** approx. adjoint and show the result */
  NFFT(adjoint)(&p);
  for (int i = 0; i < 7; i++)
    Rprintf("Here is nfft vector f_hat: %f + %f i\n", crealf(p.f_hat[i]),
            cimagf(p.f_hat[i]));

  /* NFFT(vpr_complex)(p.f_hat, K, "adjoint nfft, vector f_hat (first few
   * entries)"); */

  PutRNGstate();
  /** finalise the two dimensional plan */
  NFFT(finalize)(&p);
}

SEXP rnfft_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFFT(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFFT(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < p.d * p.M_total; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.d * p.M_total; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    error("'X' must be real or complex.");
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfft_precompute_one_psi(&p);

  Rcomplex *ffhat = COMPLEX(FHAT);
  for (j = 0; j < p.N_total; j++) {
    p.f_hat[j] = ffhat[j].r + I * ffhat[j].i;
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFFT(trafo)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_COMPLEX_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i].r = creal(p.f[i]);
    ret[i].i = cimag(p.f[i]);
  }
  NFFT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rndft_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFFT(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFFT(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < p.d * p.M_total; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.d * p.M_total; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    error("'X' must be real or complex.");
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfft_precompute_one_psi(&p);

  Rcomplex *ffhat = COMPLEX(FHAT);
  for (j = 0; j < p.N_total; j++) {
    p.f_hat[j] = ffhat[j].r + I * ffhat[j].i;
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFFT(trafo_direct)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_COMPLEX_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i].r = creal(p.f[i]);
    ret[i].i = cimag(p.f[i]);
  }
  NFFT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfft_adjoint_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFFT(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFFT(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < p.d * p.M_total; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.d * p.M_total; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    error("'X' must be real or complex.");
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfft_precompute_one_psi(&p);

  Rcomplex *ffhat = COMPLEX(FHAT);
  for (j = 0; j < p.M_total; j++) {
    p.f[j] = ffhat[j].r + I * ffhat[j].i;
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFFT(adjoint)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_COMPLEX_VECTOR(F, ret, p.N_total);
  for (i = 0; i < p.N_total; ++i) {
    ret[i].r = creal(p.f_hat[i]);
    ret[i].i = cimag(p.f_hat[i]);
  }
  NFFT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rndft_adjoint_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFFT(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFFT(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;
  if (CPLXSXP == TYPEOF(X)) {
    Rcomplex *xx = COMPLEX(X);
    for (i = 0; i < p.d * p.M_total; ++i) {
      p.x[i] = xx[i].r + I * xx[i].i;
    }
  } else if (REALSXP == TYPEOF(X)) {
    double *xx = REAL(X);
    for (int i = 0; i < p.d * p.M_total; ++i) {
      p.x[i] = xx[i] + I * 0;
    }
  } else {
    error("'X' must be real or complex.");
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfft_precompute_one_psi(&p);

  Rcomplex *ffhat = COMPLEX(FHAT);
  for (j = 0; j < p.M_total; j++) {
    p.f[j] = ffhat[j].r + I * ffhat[j].i;
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfft_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFFT(adjoint_direct)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_COMPLEX_VECTOR(F, ret, p.N_total);
  for (i = 0; i < p.N_total; ++i) {
    ret[i].r = creal(p.f_hat[i]);
    ret[i].i = cimag(p.f_hat[i]);
  }
  NFFT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfct_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFCT(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFCT(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;

  double *xx = REAL(X);
  for (int i = 0; i < p.d * p.M_total; ++i) {
    p.x[i] = xx[i];
  }

  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfct_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.N_total; j++) {
    p.f_hat[j] = ffhat[j];
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfct_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFCT(trafo)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_REAL_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i] = (p.f[i]);
  }
  NFCT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rndct_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFCT(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFCT(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;

  double *xx = REAL(X);
  for (int i = 0; i < p.d * p.M_total; ++i) {
    p.x[i] = xx[i];
  }

  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfct_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.N_total; j++) {
    p.f_hat[j] = ffhat[j];
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfct_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFCT(trafo_direct)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_REAL_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i] = (p.f[i]);
  }
  NFCT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfct_adjoint_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFCT(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFCT(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;

  double *xx = REAL(X);
  for (int i = 0; i < p.d * p.M_total; ++i) {
    p.x[i] = xx[i];
  }

  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfct_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.M_total; j++) {
    p.f[j] = ffhat[j];
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfct_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFCT(adjoint)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_REAL_VECTOR(F, ret, p.N_total);
  for (i = 0; i < p.N_total; ++i) {
    ret[i] = (p.f_hat[i]);
  }
  NFCT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rndct_adjoint_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFCT(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFCT(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;

  double *xx = REAL(X);
  for (int i = 0; i < p.d * p.M_total; ++i) {
    p.x[i] = xx[i];
  }

  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfct_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.M_total; j++) {
    p.f[j] = ffhat[j];
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfct_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFCT(adjoint_direct)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_REAL_VECTOR(F, ret, p.N_total);
  for (i = 0; i < p.N_total; ++i) {
    ret[i] = (p.f_hat[i]);
  }
  NFCT(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfst_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFST(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFST(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;

  double *xx = REAL(X);
  for (int i = 0; i < p.d * p.M_total; ++i) {
    p.x[i] = xx[i];
  }

  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfst_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.N_total; j++) {
    p.f_hat[j] = ffhat[j];
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfst_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFST(trafo)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_REAL_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i] = (p.f[i]);
  }
  NFST(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rndst_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFST(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFST(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;

  double *xx = REAL(X);
  for (int i = 0; i < p.d * p.M_total; ++i) {
    p.x[i] = xx[i];
  }

  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfst_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.N_total; j++) {
    p.f_hat[j] = ffhat[j];
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfst_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFST(trafo_direct)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_REAL_VECTOR(F, ret, p.M_total);
  for (i = 0; i < p.M_total; ++i) {
    ret[i] = (p.f[i]);
  }
  NFST(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rnfst_adjoint_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFST(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFST(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;

  double *xx = REAL(X);
  for (int i = 0; i < p.d * p.M_total; ++i) {
    p.x[i] = xx[i];
  }

  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfst_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.M_total; j++) {
    p.f[j] = ffhat[j];
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfst_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFST(adjoint)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_REAL_VECTOR(F, ret, p.N_total);
  for (i = 0; i < p.N_total; ++i) {
    ret[i] = (p.f_hat[i]);
  }
  NFST(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}

SEXP rndst_adjoint_2d(SEXP X, SEXP FHAT, SEXP M, SEXP N0, SEXP N1) {

  NFST(plan) p;
  int j;
  int m = asInteger(M);
  int n0 = asInteger(N0);
  int n1 = asInteger(N1);
  const char *error_str;
  NFST(init_2d)(&p, n0, n1, m);
  /* Rprintf("Here is N.total: %d\n", p.N_total); */
  int i = 0;

  double *xx = REAL(X);
  for (int i = 0; i < p.d * p.M_total; ++i) {
    p.x[i] = xx[i];
  }

  /* for(int i = 0; i < 7; i++) Rprintf("Here are the x[i]: %f\n",p.x[i]); */
  nfst_precompute_one_psi(&p);

  double *ffhat = REAL(FHAT);
  for (j = 0; j < p.M_total; j++) {
    p.f[j] = ffhat[j];
  }
  /* for(int i = 0; i < 7; i++) Rprintf("Here is f_hat: %f + %f
   * i\n",crealf(p.f_hat[i]), cimagf(p.f_hat[i])); */

  error_str = nfst_check(&p);
  if (error_str != 0) {
    Rprintf("Error in nfft module,\n");
    return R_NilValue;
  }

  /** trafo and show the result */
  NFST(adjoint_direct)(&p);

  /* for(int i = 0; i < 7; i++) Rprintf("Here is nfft %f + %f
   * i\n",crealf(p.f[i]),cimagf(p.f[i])); */
  /* Rprintf("Here is M_total: %d\n", p.M_total); */
  ALLOC_REAL_VECTOR(F, ret, p.N_total);
  for (i = 0; i < p.N_total; ++i) {
    ret[i] = (p.f_hat[i]);
  }
  NFST(finalize)(&p);

  UNPROTECT(1); /* s_ret */
  return F;
}
