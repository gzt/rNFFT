
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rNFFT

<!-- badges: start -->

[![R build
status](https://github.com/gzt/rNFFT/workflows/R-CMD-check/badge.svg)](https://github.com/gzt/rNFFT/actions)
[![Build
Status](https://travis-ci.org/gzt/rNFFT.svg?branch=master)](https://travis-ci.org/gzt/rNFFT)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

Writing an R wrapper for the NFFT library (nonequispaced nodes Fast
Fourier Transform)

See: \* <https://github.com/NFFT/nfft> \*
<https://tu-chemnitz.de/~potts/nfft/>

Currently, this is not in anything approaching a working state. I think
the issue may be that I need to handle this with external pointers (as
`fftw` does but `fftwtools` doesn’t).

NEWS (07/21/19): I fixed a thing so now the test function works. From
here, completing the package should be straightforward once I put a
little work in (or if anybody looks at the repo, they can take over from
here).

NEWS (07/22/19): I have implemented the 1D functions, including the
cosine and sine transforms. I don’t have the 2D functions or the 3D
functions. I do not have the ability to alter the flags in there yet,
either. However, they can be specified in the future. Caveat emptor.
EDIT: I have this working on Travis now, the issue was telling it to run
`autoconf` (I thought the R build process did that automatically, but it
doesn’t).
