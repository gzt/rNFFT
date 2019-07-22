
[![Build Status](https://travis-ci.org/gzt/rNFFT.svg?branch=master)](https://travis-ci.org/gzt/rNFFT)

# rNFFT
Writing an R wrapper for the NFFT library (nonequispaced nodes Fast Fourier Transform)

See: 
* https://github.com/NFFT/nfft
* https://tu-chemnitz.de/~potts/nfft/

Currently, this is not in anything approaching a working state. I think the issue may be 
that I need to handle this with external pointers (as `fftw` does but `fftwtools` doesn't).

NEWS (07/21/19): I fixed a thing so now the test function works. From here,
completing the package should be straightforward once I put a little work in
(or if anybody looks at the repo, they can take over from here). 

NEWS (07/22/19): I have implemented the 1D functions. I don't have the solver, nor the 2D 
functions, nor do I have this working on Travis, but it works on both my (Fedora) 
machines. I do not have the ability to alter the flags in there yet, either.
However, they can be specified in the future. Caveat emptor. EDIT: I have this working
on Travis now, the issue was telling it to run `autoconf` (I thought the R build 
process did that automatically, but it doesn't).
