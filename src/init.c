#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP test_function(SEXP, SEXP);
extern SEXP rndft_1d(SEXP , SEXP  , SEXP , SEXP );
extern SEXP rnfft_1d(SEXP , SEXP  , SEXP , SEXP );
extern SEXP rnfft_adjoint_1d(SEXP , SEXP  , SEXP , SEXP );
extern SEXP rndft_adjoint_1d(SEXP , SEXP  , SEXP , SEXP );
extern SEXP solvetest(SEXP, SEXP, SEXP);
extern SEXP rnfft_solver_1d(SEXP , SEXP , SEXP , SEXP , SEXP, SEXP );

static const R_CallMethodDef CallEntries[] = {
    {"test_function", (DL_FUNC) &test_function, 2},
    {"rndft_1d", (DL_FUNC) &rndft_1d, 4},
    {"rnfft_1d", (DL_FUNC) &rnfft_1d, 4},
    {"rnfft_adjoint_1d", (DL_FUNC) &rnfft_adjoint_1d, 4},
    {"rndft_adjoint_1d", (DL_FUNC) &rndft_adjoint_1d, 4},
    {"solvetest", (DL_FUNC) &solvetest, 3},
     {"rnfft_solver_1d", (DL_FUNC) &rnfft_solver_1d, 6},
    {NULL, NULL, 0}
};

void R_init_rNFFT(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
