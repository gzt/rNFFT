
##' 1-D Non-Uniform Fast Cosine/Sine Tranform
##'
##' Similar to the NFFT, except the cosine/sine transform instead of the
##' FFT. Note that all values of \code{x} must be in \eqn{[0,1/2)}.
##' @inheritParams ndft_1d
##' @export
##' @examples
##'
##' set.seed(20190722)
##' x <- runif(19) * 0.5
##' f_hat = runif(14)
##' 
##' nfct_1d(x, f_hat)
##' nfst_1d(x, f_hat)
##' f <- nfct_1d(x, f_hat)
##' g <- nfst_1d(x, f_hat)
##' nfct_adjoint_1d(x,f,14)
##' nfst_adjoint_1d(x,g,14)
nfct_1d <- function(x, f_hat){
    if(!(all(x <= 0.5))) stop("x must be less than 0.5")
    if(!(all(x >=0))) stop("x must be greater than 0.")
       
    M = length(x)
    N = length(f_hat)
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    #if(M < N) stop("length of X must be greater than length of f_hat")
    
    f = .Call("rnfct_1d", x, (f_hat), as.integer(M), as.integer(N), PACKAGE="rNFFT")
    return(f)
}


##' 1-D Non-Uniform Fast Cosine Tranform (Adjoint)
##' @param f frequencies for adjoint, same length as \code{x}
##' @param N number of frequencies for transform, specified for adjoint.
##' @describeIn nfct_1d The adjoint transform
##' 
##' @export
nfct_adjoint_1d <- function(x, f, N){
    M = length(x)
    if(length(f) != length(x)) stop("f must have the same length as x")
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    #if(M < N) stop("length of X must be greater than length of f_hat")
    
    fhat = .Call("rnfct_adjoint_1d", x, f, as.integer(M), as.integer(N), PACKAGE="rNFFT")
    return(fhat)
}





##' @describeIn nfct_1d The sine transform
##' @export
nfst_1d <- function(x, f_hat){
    if(!(all(x <= 0.5))) stop("x must be less than 0.5")
    if(!(all(x >=0))) stop("x must be greater than 0.")
       
    M = length(x)
    N = length(f_hat)
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    #if(M < N) stop("length of X must be greater than length of f_hat")
    
    f = .Call("rnfst_1d", x, (f_hat), as.integer(M), as.integer(N), PACKAGE="rNFFT")
    return(f)
}


##' @describeIn nfct_1d The adjoint sine transform
##' @export
nfst_adjoint_1d <- function(x, f, N){
    M = length(x)
    if(length(f) != length(x)) stop("f must have the same length as x")
    if(N%%2 != 0) stop("Must have an even number of frequencies")
    
    #if(M < N) stop("length of X must be greater than length of f_hat")
    
    fhat = .Call("rnfst_adjoint_1d", x, f, as.integer(M), as.integer(N), PACKAGE="rNFFT")
    return(fhat)
}


