
##' 2-D Non-equispaced Fourier Tranform
##'
##' The non-equispaced Fourier transform takes non-uniform samples \eqn{x}
##' from the $d$-dimensional torus \eqn{[0.5,0.5)^d}.
##'
##' The NDFT functions compute the Fourier transform directly. This is slow.
##' The NFFT functions use the FFT to compute this, which should be faster.
##' The adjoint, in this case, is not the same as the inverse. Solving the
##' inverse problem requires approximations. Here we present the 2D NDFT,
##' NFFT, and their adjoints. You most likely want to use the \code{nfft_2d}
##' and \code{nfft_adjoint_2d} functions rather than the \code{dft} functions.
##'
##' @title 2-D NFFT
##' @export
##' @param x two dimensional complex vector in \eqn{[-0.5,0.5)^2}
##' @param f_hat set of frequencies
##' @examples 
##' set.seed(20190728)
##' x <- matrix(runif(2*32*14)-.5, ncol=2, byrow = TRUE)
##' 
##' f_hat = 1:(32*14)
##' for(i in 1:(32*14)) f_hat[i] = runif(1)*1i + runif(1)
##' f_hatmatrix = matrix(f_hat, nrow = 32) 
##' f_vector <- nfft_2d(x, f_hatmatrix)
##' fd_vector <- ndft_2d(x, f_hatmatrix)
##' f_vector[1:3]
##' fd_vector[1:3]
##' newf_hat <- nfft_adjoint_2d(x,f_vector,32,14)
##' newf_hat[1,1:7]
##' newdf_hat <- ndft_adjoint_2d(x,f_vector,32,14)
##' newdf_hat[1,1:7]
nfft_2d <- function(x, f_hat){
    dims = dim(x)
    fdims = dim(f_hat)
    if(length(dims)!=2) stop("x must be 2-dimensional")
    if(length(fdims)!=2) stop("f_hat must be 2-dimensional")
    M = dims[1]
    N0 = fdims[1]
    N1 = fdims[2]
    if((N0%%2 != 0) || N1%%2 !=0) stop("Must have an even number of frequencies")
    xvec = c(t(x))
    f_hatvec = c(f_hat)
    .Call("rnfft_2d", xvec, as.complex(f_hatvec), as.integer(M), as.integer(N0), as.integer(N1), PACKAGE="rNFFT")
}

##' 2-D direct Fourier Transform
##' @export
##' @describeIn nfft_2d 
ndft_2d <- function(x, f_hat){
    dims = dim(x)
    fdims = dim(f_hat)
    if(length(dims)!=2) stop("x must be 2-dimensional")
    if(length(fdims)!=2) stop("f_hat must be 2-dimensional")
    M = dims[1]
    N0 = fdims[1]
    N1 = fdims[2]
    if((N0%%2 != 0) || N1%%2 !=0) stop("Must have an even number of frequencies")
    xvec = c(t(x))
    f_hatvec = c(f_hat)
    .Call("rndft_2d", xvec, as.complex(f_hatvec), as.integer(M), as.integer(N0), as.integer(N1), PACKAGE="rNFFT")
}



##' 2-D Non-Uniform Fast Fourier Tranform (Adjoint)
##' @param f frequencies for adjoint, same length as \code{x}
##' @param N0 number of frequencies in the first dimension for transform, specified for adjoint.
##' @param N1 number of frequencies in the second dimension for transform, specified for adjoint.
##' @describeIn nfft_2d
##' 
##' @export
nfft_adjoint_2d <- function(x, f, N0, N1){
    dims = dim(x)
    fdims = length(f)

    M = dims[1]
    if(length(dims)!=2) stop("x must be 2-dimensional")
    if(fdims != M) stop("f must have same length as x")
    
    if((N0%%2 != 0) || N1%%2 !=0) stop("Must have an even number of frequencies")

    #if(fdims[1] !=dims[1]) stop("f must have the same length as x")
      
    #if(M < N) stop("length of X must be greater than length of f_hat")
    
    fhat = .Call("rnfft_adjoint_2d", c(t(x)), as.complex(f), as.integer(M), as.integer(N0),as.integer(N1), PACKAGE="rNFFT")
    return(matrix(fhat, nrow=N0, byrow=TRUE))
}


##' 2-D Non-Uniform Direct Fourier Tranform (Adjoint)
##' @describeIn nfft_2d
##' 
##' @export
ndft_adjoint_2d <- function(x, f, N0, N1){
    dims = dim(x)
    fdims = length(f)

    M = dims[1]
    if(length(dims)!=2) stop("x must be 2-dimensional")
    if(fdims != M) stop("f must have same length as x")
    
    if((N0%%2 != 0) || N1%%2 !=0) stop("Must have an even number of frequencies")

    #if(fdims[1] !=dims[1]) stop("f must have the same length as x")
      
    #if(M < N) stop("length of X must be greater than length of f_hat")
    
    fhat = .Call("rndft_adjoint_2d", c(t(x)), as.complex(f), as.integer(M), as.integer(N0),as.integer(N1), PACKAGE="rNFFT")
    return(matrix(fhat, nrow=N0, byrow=TRUE))
}
