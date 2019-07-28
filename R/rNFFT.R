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
##' 
##' @return TRUE
##' @export
##' @author Geoffrey Z. Thompson
##'
##' @examples
##' set.seed(20190722)
##' test_nfft(19L,14L)
##' test_solve(8L,16L,9L)
test_nfft<- function(M, N){
    
  .Call("test_function", as.integer(M), as.integer(N), PACKAGE="rNFFT")

  return(TRUE)    
}

##' @export
##' @describeIn test_nfft Test solver: displays the 1D solver.
##' @param iterations Number of iterations (for solver)
test_solve <- function(M, N, iterations){
    .Call("solvetest", as.integer(M), as.integer(N), as.integer(iterations), PACKAGE = "rNFFT")
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
##'     (of length \code{N}) and length must be even.
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
    
    #if(M < N) stop("length of X must be greater than length of f_hat")
    
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
    
    #if(M < N) stop("length of X must be greater than length of f_hat")
    
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







##' 1-D NFFT solver
##'
##' Unlike in the case of the equispaced Fourier transform, the adjoint is NOT
##' the inverse of the non-equispaced Fourier transform. The solution must be found
##' numerically. 
##' @title Inverse of the 1-D NFFT.
##' @param x input vector (real, length \code{M})
##' @param f Outputted Fourier coefficients (complex, length \code{M})
##' @param N Number of coefficients for the inverse - must be even.
##' @param eps convergence criterion
##' @param iterations number of iterations for solve to run (at most)
##' @return vector of length \code{N} of solutions \eqn{\hat{f}}{f_hat}.
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
nfft_solver_1d <- function(x, f, N, eps = 1e-12, iterations = 10){
    M = length(x)
    if(length(f) != length(x)) stop("f must have the same length as x")
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    #if(M > N) stop("length of f_hat must be greater than length of X")
    
    fhat = .Call("rnfft_solver_1d", x, as.complex(f), as.integer(M), as.integer(N), eps, as.integer(iterations), PACKAGE="rNFFT")
    return(fhat)
    
}

#' Radon transform using NFFT
#'
#' This doesn't work very well. 
#' @export
#' @param image square image
#' @param Theta Number of theta slices
#' @param Rho Number of rho slices
#' @param fn Whether to use polar or linotype (polar by default)
#' @examples
#' P <- PET::phantom()
#'
#' P_radon <- nfft_radon(P, 514, 514, fn = "polar")
#' #image(P_radon)
#' P_inv <- nfft_inv_radon((P_radon), N = 257,iter = 5, fn = "polar")
#' #image(P_inv)
nfft_radon <- function(image, Theta = 181, Rho = 2*round(sqrt(sum((dim(image))^2))/2)+1, fn = "polar"){
   ## image must be N x N or a vector NxN
    dims = dim(image)
    if(length(dims) != 2) stop("Not a 2D image")
    
    N = dims[1]
  

    if(dims[1] != dims[2]) stop("Image not square.")

    fntag = ifelse(fn == "polar", 1, 0)
    

    ret_matrix = .Call("c_radon", c(image), as.integer(fntag), as.integer(N), as.integer(Theta), as.integer(Rho))

    matrix(ret_matrix, Rho, Theta, byrow = FALSE)
    
}

#' @describeIn nfft_radon Inverse Radon Transform using NFFT
#' @param N size of image
#' @param iter number of iterations for inverse
#' @export
nfft_inv_radon <- function(image, N, iter = 10, fn = "polar"){
   ## image must be N x N or a vector NxN
    dims = dim(image)
    if(length(dims) != 2) stop("Not a 2D image")
    
    Rho = dims[1]
    Theta = dims[2]

    fntag = ifelse(fn == "polar", 1, 0)
    

    ret_matrix = .Call("c_inv_radon", c((image)), as.integer(fntag), as.integer(N), as.integer(Theta), as.integer(Rho), as.integer(iter))

    matrix(ret_matrix, N)
    
}


.onUnload <- function(libpath) {
 library.dynam.unload("rNFFT", libpath)
}
