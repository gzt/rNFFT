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
##' Generates Test Output. Compare the output the results
##'
##' @title Test Function
##' @param M integer, number of nodes
##' @param N integer, number of frequences (even, less than \code{M}).
##' 
##' @return TRUE
##' @export
##' @author Geoffrey Z. Thompson
##'
##' @examples
##' set.seed(20190722)
##' test_nfft(19L,14L)
##' 
test_nfft<- function(M, N){
    
  .Call("test_function", as.integer(M), as.integer(N), PACKAGE="rNFFT")

  return(TRUE)    
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
##' NFFT, and their adjoints. You most likely want to use the \code{nfft_1d}
##' and \code{nfft_adjoint_1d} functions rather than the \code{dft} functions.
##'
##' 
##' @title 1-D NFFT
##' @param x (real) vector of nodes of length \code{M}, must be in \code{[-0.5,0.5)}.
##' @param f_hat (complex) vector of \eqn{\hat{f}}{f_hat} entries -
##'     must be shorter than \code{x} (of length \code{N}) and length must be even.
##' @return vector of \eqn{f}, the results of the transform (of length \code{M}).
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
ndft_1d <- function(x, f_hat){
    M = length(x)
    N = length(f_hat)
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    if(M < N) stop("length of X must be greater than length of f_hat")
    
    f = .Call("rndft_1d", x, as.complex(f_hat), as.integer(M), as.integer(N), PACKAGE="rNFFT")
    return(f)
}





##' 1-D Non-Uniform Fast Fourier Tranform
##' @describeIn ndft_1d
##' 
##' @export
nfft_1d <- function(x, f_hat){
    M = length(x)
    N = length(f_hat)
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    if(M < N) stop("length of X must be greater than length of f_hat")
    
    f = .Call("rnfft_1d", x, as.complex(f_hat), as.integer(M), as.integer(N), PACKAGE="rNFFT")
    return(f)
}


##' 1-D Non-Uniform Fast Fourier Tranform (Adjoint)
##' @param f frequencies for adjoint, same length as \code{x}
##' @param N number of frequencies for transform, specified for adjoint.
##' @describeIn ndft_1d
##' 
##' @export
nfft_adjoint_1d <- function(x, f, N){
    M = length(x)
    if(length(f) != length(x)) stop("f must have the same length as x")
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    if(M < N) stop("length of X must be greater than length of f_hat")
    
    fhat = .Call("rnfft_adjoint_1d", x, as.complex(f), as.integer(M), as.integer(N), PACKAGE="rNFFT")
    return(fhat)
}


##' 1-D Non-Uniform Direct Fourier Tranform
##' @describeIn ndft_1d
##' 
##' @export
ndft_adjoint_1d <- function(x, f, N){
    M = length(x)
    if(length(f) != length(x)) stop("f must have the same length as x")
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    if(M < N) stop("length of X must be greater than length of f_hat")
    
    fhat = .Call("rndft_adjoint_1d", x, as.complex(f), as.integer(M), as.integer(N), PACKAGE="rNFFT")
    return(fhat)
}

.onUnload <- function(libpath) {
 library.dynam.unload("rNFFT", libpath)
}


