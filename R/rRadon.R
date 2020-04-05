
#' Radon transform using NFFT
#'
#' A radon transformation using the non-equispaced
#' fast Fourier transform (NFFT).
#' Requires a square image.
#' @export
#' @param image square image
#' @param theta Number of theta slices
#' @param rho Number of rho slices
#' @param fn Whether to use polar or linotype (polar by default)
#' @examples
#' P <- PET::phantom()
#' image(P)
#' P_radon <- nfft_radon(P, 514, 514, fn = "polar")
#' image(P_radon)
#' P_inv <- nfft_inv_radon((P_radon), n = 257, iter = 2, fn = "polar")
#' image(P_inv)
nfft_radon <- function(image,
                       theta = 181,
                       rho = 2 * round(sqrt(sum((dim(image))^2)) / 2) + 1,
                       fn = "polar") {
  ## image must be n x n or a vector nxn
  dims <- dim(image)
  if (length(dims) != 2) stop("Not a 2D image")

  n <- dims[1]


  if (dims[1] != dims[2]) stop("Image not square.")

  fntag <- ifelse(fn == "polar", 1, 0)


  ret_matrix <- .Call(
    "c_radon", c(image), as.integer(fntag), as.integer(n),
    as.integer(theta), as.integer(rho)
  )

  matrix(ret_matrix, rho, theta, byrow = FALSE)
}

#' @describeIn nfft_radon Inverse Radon Transform using NFFT
#' @param n size of image
#' @param iter number of iterations for inverse
#' @export
nfft_inv_radon <- function(image, n, iter = 10, fn = "polar") {
  ## image must be n x n or a vector nxn
  dims <- dim(image)
  if (length(dims) != 2) stop("Not a 2D image")

  rho <- dims[1]
  theta <- dims[2]

  fntag <- ifelse(fn == "polar", 1, 0)


  ret_matrix <- .Call(
    "c_inv_radon", c((image)), as.integer(fntag),
    as.integer(n), as.integer(theta), as.integer(rho),
    as.integer(iter)
  )

  matrix(ret_matrix, n)
}
