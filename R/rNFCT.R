
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
##' nfst_adjoint_1d(x,g,14) # note the length!
nfct_1d <- function(x, f_hat) {
  if (!(all(x <= 0.5))) stop("x must be less than 0.5")
  if (!(all(x >= 0))) stop("x must be greater than 0.")

  m <- length(x)
  n <- length(f_hat)
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  f <- .Call("rnfct_1d", x, (f_hat), as.integer(m), as.integer(n),
    PACKAGE = "rNFFT"
  )
  return(f)
}


##' 1-D Non-Uniform Fast Cosine Tranform (Adjoint)
##' @param f frequencies for adjoint, same length as \code{x}
##' @param n number of frequencies for transform, specified for adjoint.
##' @describeIn nfct_1d The adjoint transform
##'
##' @export
nfct_adjoint_1d <- function(x, f, n) {
  m <- length(x)
  if (length(f) != length(x)) stop("f must have the same length as x")
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  fhat <- .Call("rnfct_adjoint_1d", x, f, as.integer(m), as.integer(n),
    PACKAGE = "rNFFT"
  )
  return(fhat)
}





##' @describeIn nfct_1d The sine transform
##' @export
nfst_1d <- function(x, f_hat) {
  if (!(all(x <= 0.5))) stop("x must be less than 0.5")
  if (!(all(x >= 0))) stop("x must be greater than 0.")

  m <- length(x)
  n <- length(f_hat)
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  f <- .Call("rnfst_1d", x, (f_hat), as.integer(m), as.integer(n),
    PACKAGE = "rNFFT"
  )
  return(f)
}


##' @describeIn nfct_1d The adjoint sine transform
##' @export
nfst_adjoint_1d <- function(x, f, n) {
  m <- length(x)
  if (length(f) != length(x)) stop("f must have the same length as x")
  if (n %% 2 != 0) stop("Must have an even number of frequencies")

  fhat <- .Call("rnfst_adjoint_1d", x, f, as.integer(m), as.integer(n),
    PACKAGE = "rNFFT"
  )
  return(fhat)
}



##' 2-D Non-equispaced Trigonometric Tranform
##'
##' @description The non-equispaced trigonometric transforms
##' (sine and cosine) take non-uniform samples \eqn{x} from the $d$-dimensional
##' torus \eqn{[0,0.5)^d}.
##'
##' The NDCT functions compute the cosine transform directly. This is slow.
##' The NFCT functions use the FFT to compute this, which should be faster.
##' The adjoint, in this case, is not the same as the inverse. Solving the
##' inverse problem requires approximations. Here we present the 2D NDCT,
##' NFCT, and their adjoints. You most likely want to use the \code{nfct_2d}
##' and \code{nfct_adjoint_2d} functions rather than the \code{dct} functions.
##'
##' @title 2-D NFCT/NFST
##' @export
##' @param x two dimensional complex vector in \eqn{[0,0.5)^2}
##' @param f_hat set of frequencies
##' @examples
##' set.seed(20190728)
##' x <- matrix(runif(2*32*14)*.5, ncol=2, byrow = TRUE)
##'
##' f_hat = 1:(32*14)
##' for(i in 1:(32*14)) f_hat[i] = runif(1)
##' f_hatmatrix = matrix(f_hat, nrow = 32)
##' f_vector <- nfct_2d(x, f_hatmatrix)
##' fs_vector <-  nfst_2d(x, f_hatmatrix)
##' fd_vector <- ndct_2d(x, f_hatmatrix)
##' fsd_vector <-  ndst_2d(x, f_hatmatrix)
##' f_vector[1:3]
##' fd_vector[1:3]
##' fs_vector[1:3]
##' fsd_vector[1:3]
##' newf_hat <- nfct_adjoint_2d(x,f_vector,32,14)
##' newf_hat[1,1:7]
##' newdf_hat <- ndct_adjoint_2d(x,f_vector,32,14)
##' newdf_hat[1,1:7]
##' newfs_hat <- nfst_adjoint_2d(x,f_vector,32,14)
##' newfs_hat[1,1:7]
##' newdsf_hat <- ndst_adjoint_2d(x,f_vector,32,14)
##' newdsf_hat[1,1:7]
nfct_2d <- function(x, f_hat) {
  dims <- dim(x)
  fdims <- dim(f_hat)
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (length(fdims) != 2) stop("f_hat must be 2-dimensional")
  m <- dims[1]
  n0 <- fdims[1]
  n1 <- fdims[2]
  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }
  xvec <- c(t(x))
  f_hatvec <- c(f_hat)
  .Call("rnfct_2d", xvec, (f_hatvec), as.integer(m), as.integer(n0),
    as.integer(n1),
    PACKAGE = "rNFFT"
  )
}

##' 2-D direct cosine Transform
##'
##' @describeIn nfct_2d
##'
##' @export
ndct_2d <- function(x, f_hat) {
  dims <- dim(x)
  fdims <- dim(f_hat)
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (length(fdims) != 2) stop("f_hat must be 2-dimensional")
  m <- dims[1]
  n0 <- fdims[1]
  n1 <- fdims[2]
  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }
  xvec <- c(t(x))
  f_hatvec <- c(f_hat)
  .Call("rndct_2d", xvec, (f_hatvec), as.integer(m), as.integer(n0),
    as.integer(n1),
    PACKAGE = "rNFFT"
  )
}



##' 2-D Non-Uniform Fast Cosine Tranform (Adjoint)
##' @param f frequencies for adjoint, same length as \code{x}
##' @param n0 number of frequencies in the first dimension for transform,
##'     specified for adjoint.
##' @param n1 number of frequencies in the second dimension for transform,
##'     specified for adjoint.
##' @describeIn nfct_2d
##'
##' @export
nfct_adjoint_2d <- function(x, f, n0, n1) {
  dims <- dim(x)
  fdims <- length(f)

  m <- dims[1]
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (fdims != m) stop("f must have same length as x")

  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }

  fhat <- .Call("rnfct_adjoint_2d", c(t(x)), (f), as.integer(m),
    as.integer(n0), as.integer(n1),
    PACKAGE = "rNFFT"
  )
  return(matrix(fhat, nrow = n0, byrow = TRUE))
}


##' 2-D Non-Uniform Direct Cosine Tranform (Adjoint)
##' @describeIn nfct_2d
##'
##' @export
ndct_adjoint_2d <- function(x, f, n0, n1) {
  dims <- dim(x)
  fdims <- length(f)

  m <- dims[1]
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (fdims != m) stop("f must have same length as x")

  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }

  fhat <- .Call("rndct_adjoint_2d", c(t(x)), (f), as.integer(m),
    as.integer(n0), as.integer(n1),
    PACKAGE = "rNFFT"
  )
  return(matrix(fhat, nrow = n0, byrow = TRUE))
}


##' 2-D Non-Uniform Fast Sine Transform
##' @describeIn nfct_2d
##'
##' @export
nfst_2d <- function(x, f_hat) {
  dims <- dim(x)
  fdims <- dim(f_hat)
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (length(fdims) != 2) stop("f_hat must be 2-dimensional")
  m <- dims[1]
  n0 <- fdims[1]
  n1 <- fdims[2]
  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }
  xvec <- c(t(x))
  f_hatvec <- c(f_hat)
  .Call("rnfst_2d", xvec, (f_hatvec), as.integer(m), as.integer(n0),
    as.integer(n1),
    PACKAGE = "rNFFT"
  )
}

##' 2-D direct sine Transform
##' @describeIn nfct_2d
##'
##' @export
ndst_2d <- function(x, f_hat) {
  dims <- dim(x)
  fdims <- dim(f_hat)
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (length(fdims) != 2) stop("f_hat must be 2-dimensional")
  m <- dims[1]
  n0 <- fdims[1]
  n1 <- fdims[2]
  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }
  xvec <- c(t(x))
  f_hatvec <- c(f_hat)
  .Call("rndst_2d", xvec, (f_hatvec), as.integer(m), as.integer(n0),
    as.integer(n1),
    PACKAGE = "rNFFT"
  )
}



##' 2-D Non-Uniform Fast Sine Tranform (Adjoint)
##' @describeIn nfct_2d
##'
##' @export
nfst_adjoint_2d <- function(x, f, n0, n1) {
  dims <- dim(x)
  fdims <- length(f)

  m <- dims[1]
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (fdims != m) stop("f must have same length as x")

  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }

  fhat <- .Call("rnfst_adjoint_2d", c(t(x)), (f), as.integer(m),
    as.integer(n0), as.integer(n1),
    PACKAGE = "rNFFT"
  )
  return(matrix(fhat, nrow = (n0 - 1), byrow = TRUE))
}


##' 2-D Non-Uniform Direct Sine Tranform (Adjoint)
##' @describeIn nfct_2d
##'
##' @export
ndst_adjoint_2d <- function(x, f, n0, n1) {
  dims <- dim(x)
  fdims <- length(f)

  m <- dims[1]
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (fdims != m) stop("f must have same length as x")

  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }


  fhat <- .Call("rndst_adjoint_2d", c(t(x)), (f), as.integer(m),
    as.integer(n0), as.integer(n1),
    PACKAGE = "rNFFT"
  )
  return(matrix(fhat, nrow = (n0 - 1), byrow = TRUE))
}
