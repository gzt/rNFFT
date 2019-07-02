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
##' Generates Test Output
##' @title Test Function
##' @param M integer, number of nodes
##' @param N integer
##' @return Null
##' @export
##' @author Geoffrey Z. Thompson
##'
##' @examples
##' 
##' test(19L,14L)
##' 
test <- function(M, N){
  .Call("test_function", as.integer(M), as.integer(N), PACKAGE="rNFFT")

    
}


.onUnload <- function(libpath) {
 library.dynam.unload("rNFFT", libpath)
}
