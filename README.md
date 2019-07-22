
[![Build Status](https://travis-ci.org/gzt/rNFFT.svg?branch=master)](https://travis-ci.org/gzt/rNFFT)

# rNFFT
Writing an R wrapper for the NFFT library (nonequispaced nodes Fast Fourier Transform)

See: 
* https://github.com/NFFT/nfft
* https://tu-chemnitz.de/~potts/nfft/

Currently, this is not in anything approaching a working state. I think the issue may be that I need to handle this with external pointers (as `fftw` does but `fftwtools` doesn't).

I may give up on this and set my sights on [FINUFFT](https://finufft.readthedocs.io/en/latest/index.html) with Rcpp, since I am having no luck in getting this to work. It's puzzling.

NEWS (07/21/19): I fixed a thing so now the test function works. From here,
completing the package should be straightforward once I put a little work in
(or if anybody looks at the repo, they can take over from here). 
