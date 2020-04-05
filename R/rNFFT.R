##
## rNFFT.R - glue between user and rnfft.c
##
## Authors:
##  Geoffrey Thompson <gzthompson@gmail.com>


#' R wrapper for the NFFT library (nonequispaced nodes Fast Fourier Transform)
#'
#' Provides a wrapper around the NFFT 3.0 C library, which implements
#' the fast Fourier transform for non-equispaced nodes, its adjoint, and its
#' inversion. This provides functions for the NFFT and Solver routines.
#'
#' @docType package
#' @name rNFFT
#' @useDynLib rNFFT, .registration = TRUE
NULL

##' Test Function
##'
##' Generates Test Output specified in NFFT documentation and source code.
##' Compare the output to the results of the examples.
##'
##' @title Test Function
##' @param m integer, number of nodes
##' @param n integer, number of frequences (even, less than `m`).
##'
##'
##' @return TRUE
##' @export
##' @author Geoffrey Z. Thompson
##'
##' @examples
##' set.seed(20190722)
##' test_nfft(19L,14L)
##' test_solve(8L,16L,9L)
##' set.seed(20190728)
##' nfft_test_2d()
test_nfft <- function(m, n) {
  .Call("test_function", as.integer(m), as.integer(n), PACKAGE = "rNFFT")

  return(TRUE)
}

##' @export
##' @describeIn test_nfft Test solver: displays the 1D solver.
##' @param iterations Number of iterations (for solver)
test_solve <- function(m, n, iterations) {
  .Call("solvetest", as.integer(m), as.integer(n), as.integer(iterations),
    PACKAGE = "rNFFT"
  )
  return(TRUE)
}

##' @export
##' @describeIn test_nfft 2D test: displays the 2D test (no input)
nfft_test_2d <- function() {
  .Call("nfft_2dtest", PACKAGE = "rNFFT")
  TRUE
}

##' 1-D Non-Uniform Direct Fourier Tranform
##'
##' The non-uniform Fourier transform takes non-uniform samples \eqn{x}
##' from the $d$-dimensional torus \eqn{[0.5,0.5)^d}.
##'
##' The NDFT functions compute the Fourier transform directly. This is slow.
##' The NFFT functions use the FFT to compute this, which should be faster.
##' The adjoint, in this case, is not the same as the inverse. Solving the
##' inverse problem requires approximations. Here we present the 1D NDFT,
##' NFFT, and their adjoints. You most likely want to use the `nfft_1d`
##' and `nfft_adjoint_1d` functions rather than the `dft` functions.
##'
##'
##' @title 1-D NFFT
##' @param x (real) vector of nodes of length `m`,
##'    must be in `[-0.5,0.5)`.
##' @param f_hat (complex) vector of \eqn{\hat{f}}{f_hat} entries -
##'     (of length `n`) and length must be even.
##' @return vector of \eqn{f}, the results of the transform
##'     (of length `m`).
##'
##' @references Keiner, J., Kunis, S., and Potts, D. ''Using NFFT 3 - a software
##'     library for various nonequispaced fast Fourier transforms'' ACM Trans.
##'     Math. Software,36, Article 19, 1-30, 2009.
##' @author Geoffrey Z. Thompson
##' @export
##' @examples
##' set.seed(20190722)
##' x <- runif(19) - 0.5
##' #f_hat <- runif(14) + 1i*runif(14)
##' f_hat = 1:14
##' for(i in 1:14) f_hat[i] = runif(1)*1i + runif(1)
##'
##' ndft_1d(x, f_hat)
##' nfft_1d(x, f_hat)
##' f <- nfft_1d(x, f_hat)
##' nfft_adjoint_1d(x,f,14)
##' ndft_adjoint_1d(x,f,14)
##'
ndft_1d <- function(x, f_hat) {
  m <- length(x)
  n <- length(f_hat)
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  if (m < n) stop("length of X must be greater than length of f_hat")

  f <- .Call("rndft_1d", x, as.complex(f_hat), as.integer(m),
    as.integer(n),
    PACKAGE = "rNFFT"
  )
  return(f)
}

##' 1-D Non-Uniform Fast Fourier Tranform
##' @describeIn ndft_1d
##'
##' @export
nfft_1d <- function(x, f_hat) {
  m <- length(x)
  n <- length(f_hat)
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  f <- .Call("rnfft_1d", x, as.complex(f_hat), as.integer(m),
    as.integer(n),
    PACKAGE = "rNFFT"
  )
  return(f)
}



##' 1-D Non-Uniform Fast Fourier Tranform (Adjoint)
##' @param f frequencies for adjoint, same length as `x`
##' @param n number of frequencies for transform, specified for adjoint.
##' @describeIn ndft_1d
##'
##' @export
nfft_adjoint_1d <- function(x, f, n) {
  m <- length(x)
  if (length(f) != length(x)) stop("f must have the same length as x")
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  fhat <- .Call("rnfft_adjoint_1d", x, as.complex(f), as.integer(m),
    as.integer(n),
    PACKAGE = "rNFFT"
  )
  return(fhat)
}





##' 1-D Non-Uniform Direct Fourier Tranform (Adjoint)
##' @describeIn ndft_1d
##'
##' @export
ndft_adjoint_1d <- function(x, f, n) {
  m <- length(x)
  if (length(f) != length(x)) stop("f must have the same length as x")
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  if (m < n) stop("length of X must be greater than length of f_hat")

  fhat <- .Call("rndft_adjoint_1d", x, as.complex(f), as.integer(m),
    as.integer(n),
    PACKAGE = "rNFFT"
  )
  return(fhat)
}

##' 1-D NFFT solver
##'
##' Unlike in the case of the equispaced Fourier transform, the adjoint is NOT
##' the inverse of the non-equispaced Fourier transform. The solution must be
##' found numerically.
##' @title Inverse of the 1-D NFFT.
##' @param x input vector (real, length `m`)
##' @param f Outputted Fourier coefficients (complex, length `m`)
##' @param n Number of coefficients for the inverse - must be even.
##' @param eps convergence criterion
##' @param iterations number of iterations for solve to run (at most)
##' @return vector of length `n` of solutions \eqn{\hat{f}}{f_hat}.
##' @export
##' @examples
##' set.seed(20190722)
##' x = runif(8) - 0.5
##' f = 1:8
##' for(i in 1:8) f[i] = runif(1)*1i + runif(1)
##' x
##' f
##' nfft_solver_1d(x,f,16,1e-35,8)
##' nfft_1d(x, nfft_solver_1d(x,f,16,eps = 1e-35))
##' @author Geoffrey Z. Thompson
nfft_solver_1d <- function(x, f, n, eps = 1e-12, iterations = 10) {
  m <- length(x)
  if (length(f) != length(x)) stop("f must have the same length as x")
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  fhat <- .Call("rnfft_solver_1d", x, as.complex(f), as.integer(m),
    as.integer(n), eps, as.integer(iterations),
    PACKAGE = "rNFFT"
  )
  return(fhat)
}

.onUnload <- function(libpath) {
  library.dynam.unload("rNFFT", libpath)
}
