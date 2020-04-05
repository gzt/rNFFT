test_that("test", {
  expect_true(test_nfft(19, 14))
  set.seed(20190722)
  x <- runif(19) - 0.5
  f_hat <- 1:14
  for (i in 1:14) f_hat[i] <- runif(1) * 1i + runif(1)
  expect_equal(ndft_1d(x, f_hat), nfft_1d(x, f_hat))
  f <- nfft_1d(x, f_hat)
  expect_equal(nfft_adjoint_1d(x, f, 14), ndft_adjoint_1d(x, f, 14))
})
