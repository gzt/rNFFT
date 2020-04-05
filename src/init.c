#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP test_function(SEXP, SEXP);
extern SEXP rndft_1d(SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfft_1d(SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfft_adjoint_1d(SEXP, SEXP, SEXP, SEXP);
extern SEXP rndft_adjoint_1d(SEXP, SEXP, SEXP, SEXP);
extern SEXP solvetest(SEXP, SEXP, SEXP);
extern SEXP rnfft_solver_1d(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfct_1d(SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfct_adjoint_1d(SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfst_1d(SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfst_adjoint_1d(SEXP, SEXP, SEXP, SEXP);
extern SEXP c_radon(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP c_inv_radon(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfft_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rndft_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern void nfft_2dtest();
extern SEXP rnfft_adjoint_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rndft_adjoint_2d(SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP rnfct_adjoint_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfst_adjoint_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rndct_adjoint_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rndst_adjoint_2d(SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP rnfct_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rnfst_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rndct_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rndst_2d(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"test_function", (DL_FUNC)&test_function, 2},
    {"rndft_1d", (DL_FUNC)&rndft_1d, 4},
    {"rnfft_1d", (DL_FUNC)&rnfft_1d, 4},
    {"rnfft_adjoint_1d", (DL_FUNC)&rnfft_adjoint_1d, 4},
    {"rndft_adjoint_1d", (DL_FUNC)&rndft_adjoint_1d, 4},
    {"solvetest", (DL_FUNC)&solvetest, 3},
    {"rnfft_solver_1d", (DL_FUNC)&rnfft_solver_1d, 6},
    {"rnfct_1d", (DL_FUNC)&rnfct_1d, 4},
    {"rnfct_adjoint_1d", (DL_FUNC)&rnfct_adjoint_1d, 4},
    {"rnfst_1d", (DL_FUNC)&rnfst_1d, 4},
    {"rnfst_adjoint_1d", (DL_FUNC)&rnfst_adjoint_1d, 4},
    {"c_radon", (DL_FUNC)&c_radon, 5},
    {"c_inv_radon", (DL_FUNC)&c_inv_radon, 6},
    {"rnfft_2d", (DL_FUNC)&rnfft_2d, 5},
    {"rndft_2d", (DL_FUNC)&rndft_2d, 5},
    {"nfft_2dtest", (DL_FUNC)&nfft_2dtest, 0},
    {"rnfft_adjoint_2d", (DL_FUNC)&rnfft_adjoint_2d, 5},
    {"rndft_adjoint_2d", (DL_FUNC)&rndft_adjoint_2d, 5},
    {"rnfct_2d", (DL_FUNC)&rnfct_2d, 5},
    {"rnfct_adjoint_2d", (DL_FUNC)&rnfct_adjoint_2d, 5},
    {"rnfst_2d", (DL_FUNC)&rnfst_2d, 5},
    {"rnfst_adjoint_2d", (DL_FUNC)&rnfst_adjoint_2d, 5},
    {"rndct_2d", (DL_FUNC)&rndct_2d, 5},
    {"rndct_adjoint_2d", (DL_FUNC)&rndct_adjoint_2d, 5},
    {"rndst_2d", (DL_FUNC)&rndst_2d, 5},
    {"rndst_adjoint_2d", (DL_FUNC)&rndst_adjoint_2d, 5},
    {NULL, NULL, 0}};

void R_init_rNFFT(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
