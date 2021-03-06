% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rNFCT.R
\name{nfct_2d}
\alias{nfct_2d}
\alias{ndct_2d}
\alias{nfct_adjoint_2d}
\alias{ndct_adjoint_2d}
\alias{nfst_2d}
\alias{ndst_2d}
\alias{nfst_adjoint_2d}
\alias{ndst_adjoint_2d}
\title{2-D NFCT/NFST}
\usage{
nfct_2d(x, f_hat)

ndct_2d(x, f_hat)

nfct_adjoint_2d(x, f, n0, n1)

ndct_adjoint_2d(x, f, n0, n1)

nfst_2d(x, f_hat)

ndst_2d(x, f_hat)

nfst_adjoint_2d(x, f, n0, n1)

ndst_adjoint_2d(x, f, n0, n1)
}
\arguments{
\item{x}{two dimensional complex vector in \eqn{[0,0.5)^2}}

\item{f_hat}{set of frequencies}

\item{f}{frequencies for adjoint, same length as \code{x}}

\item{n0}{number of frequencies in the first dimension for transform,
specified for adjoint.}

\item{n1}{number of frequencies in the second dimension for transform,
specified for adjoint.}
}
\description{
The non-equispaced trigonometric transforms
(sine and cosine) take non-uniform samples \eqn{x} from the $d$-dimensional
torus \eqn{[0,0.5)^d}.

The NDCT functions compute the cosine transform directly. This is slow.
The NFCT functions use the FFT to compute this, which should be faster.
The adjoint, in this case, is not the same as the inverse. Solving the
inverse problem requires approximations. Here we present the 2D NDCT,
NFCT, and their adjoints. You most likely want to use the \code{nfct_2d}
and \code{nfct_adjoint_2d} functions rather than the \code{dct} functions.
}
\details{
2-D Non-equispaced Trigonometric Tranform
}
\section{Functions}{
\itemize{
\item \code{ndct_2d}: 

\item \code{nfct_adjoint_2d}: 

\item \code{ndct_adjoint_2d}: 

\item \code{nfst_2d}: 

\item \code{ndst_2d}: 

\item \code{nfst_adjoint_2d}: 

\item \code{ndst_adjoint_2d}: 
}}

\examples{
set.seed(20190728)
x <- matrix(runif(2*32*14)*.5, ncol=2, byrow = TRUE)

f_hat = 1:(32*14)
for(i in 1:(32*14)) f_hat[i] = runif(1)
f_hatmatrix = matrix(f_hat, nrow = 32)
f_vector <- nfct_2d(x, f_hatmatrix)
fs_vector <-  nfst_2d(x, f_hatmatrix)
fd_vector <- ndct_2d(x, f_hatmatrix)
fsd_vector <-  ndst_2d(x, f_hatmatrix)
f_vector[1:3]
fd_vector[1:3]
fs_vector[1:3]
fsd_vector[1:3]
newf_hat <- nfct_adjoint_2d(x,f_vector,32,14)
newf_hat[1,1:7]
newdf_hat <- ndct_adjoint_2d(x,f_vector,32,14)
newdf_hat[1,1:7]
newfs_hat <- nfst_adjoint_2d(x,f_vector,32,14)
newfs_hat[1,1:7]
newdsf_hat <- ndst_adjoint_2d(x,f_vector,32,14)
newdsf_hat[1,1:7]
}
