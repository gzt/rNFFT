###   test-test-test.R
###   rNFFT : Tests for NFFT
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
