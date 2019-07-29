
#' Radon transform using NFFT
#'
#' A radon transformation using the non-equispaced fast Fourier transform (NFFT).
#' Requires a square image. 
#' @export
#' @param image square image
#' @param Theta Number of theta slices
#' @param Rho Number of rho slices
#' @param fn Whether to use polar or linotype (polar by default)
#' @examples
#' P <- PET::phantom()
#'
#' P_radon <- nfft_radon(P, 514,514, fn = "polar")
#' image(P_radon)
#' P_inv <- nfft_inv_radon((P_radon), N = 257,iter = 2, fn = "polar")
#' image(P_inv)
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

