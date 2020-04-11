###   NFFT2d.R
###   rNFFT : 2-D functions for rNFFT
###   Copyright (C) 2020  GZ Thompson <gzthompson@gmail.com>
###
###   This program is free software; you can redistribute it and/or modify
###   it under the terms of the GNU General Public License as published by
###   the Free Software Foundation; either version 3 of the License, or
###   (at your option) any later version.
###
###   This program is distributed in the hope that it will be useful,
###   but WITHOUT ANY WARRANTY; without even the implied warranty of
###   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###   GNU General Public License for more details.
###
###     You should have received a copy of the GNU General Public License
###   along with this program; if not, a copy is available at
###   https://www.R-project.org/Licenses/

##' 2-D Non-equispaced Fourier Tranform
##'
##' The non-equispaced Fourier transform takes non-uniform samples \eqn{x}
##' from the $d$-dimensional torus \eqn{[0.5,0.5)^d}.
##'
##' The NDFT functions compute the Fourier transform directly. This is slow.
##' The NFFT functions use the FFT to compute this, which should be faster.
##' The adjoint, in this case, is not the same as the inverse. Solving the
##' inverse problem requires approximations. Here we present the 2D NDFT,
##' NFFT, and their adjoints. You most likely want to use the `nfft_2d`
##' and `nfft_adjoint_2d` functions rather than the `dft` functions.
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
nfft_2d <- function(x, f_hat) {
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
  .Call("rnfft_2d", xvec, as.complex(f_hatvec), as.integer(m),
    as.integer(n0), as.integer(n1),
    PACKAGE = "rNFFT"
  )
}

##' 2-D direct Fourier Transform
##'
##' @export
##' @describeIn nfft_2d
##'
ndft_2d <- function(x, f_hat) {
  dims <- dim(x)
  fdims <- dim(f_hat)
  if (length(dims) != 2) {
    stop("x must be 2-dimensional")
  }
  if (length(fdims) != 2) {
    stop("f_hat must be 2-dimensional")
  }
  m <- dims[1]
  n0 <- fdims[1]
  n1 <- fdims[2]
  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }
  xvec <- c(t(x))
  f_hatvec <- c(f_hat)
  .Call("rndft_2d", xvec, as.complex(f_hatvec), as.integer(m),
    as.integer(n0), as.integer(n1),
    PACKAGE = "rNFFT"
  )
}



##' 2-D Non-Uniform Fast Fourier Tranform (Adjoint)
##' @param f frequencies for adjoint, same length as `x`
##' @param n0 number of frequencies in the first dimension for
##' transform, specified for adjoint.
##' @param n1 number of frequencies in the second dimension for
##' transform, specified for adjoint.
##' @describeIn nfft_2d
##'
##' @export
nfft_adjoint_2d <- function(x, f, n0, n1) {
  dims <- dim(x)
  fdims <- length(f)

  m <- dims[1]
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (fdims != m) stop("f must have same length as x")

  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }

  fhat <- .Call("rnfft_adjoint_2d", c(t(x)), as.complex(f), as.integer(m),
    as.integer(n0), as.integer(n1),
    PACKAGE = "rNFFT"
  )
  return(matrix(fhat, nrow = n0, byrow = TRUE))
}


##' 2-D Non-Uniform Direct Fourier Tranform (Adjoint)
##' @describeIn nfft_2d
##'
##' @export
ndft_adjoint_2d <- function(x, f, n0, n1) {
  dims <- dim(x)
  fdims <- length(f)

  m <- dims[1]
  if (length(dims) != 2) stop("x must be 2-dimensional")
  if (fdims != m) stop("f must have same length as x")

  if ((n0 %% 2 != 0) || n1 %% 2 != 0) {
    stop("Must have an even number of frequencies")
  }

  fhat <- .Call("rndft_adjoint_2d", c(t(x)), as.complex(f), as.integer(m),
    as.integer(n0), as.integer(n1),
    PACKAGE = "rNFFT"
  )
  return(matrix(fhat, nrow = n0, byrow = TRUE))
}
