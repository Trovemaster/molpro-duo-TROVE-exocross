# TROVE binary used for the cluster runs

* Source: Trovemaster/TROVE `develop` at commit 930f7f64635c7683c6703b5275fff62d3b7d519e.
  This commit includes the CH3OH fixes 70e668d and 1ecc400 and the `PThamiltonian_contract`
  fix 369ef18.
* Local changes:
  * `makefile`: COMPILER=gfortran, `-march=x86-64-v3`, `LAPACK = -llapack -lblas`, and
    `ddoti.f90` added to SRCS.
  * `ddoti.f90`: a portable replacement for the MKL sparse `ddoti`.
* Build environment: Ubuntu 22.04, GNU Fortran 11.4.0, libopenblas-dev 0.3.20 (pthread) and
  liblapack-dev.
* Compiler flags: `-cpp -std=gnu -fopenmp -march=x86-64-v3 -ffree-line-length-512 -fcray-pointer
  -fallow-argument-mismatch -O3`.
* Run-time libraries: liblapack.so.3 and libblas.so.3 (resolved to OpenBLAS), libgfortran.so.5,
  libgomp.so.1 and glibc ≥ 2.34.
* Validation: the P6 test basis gives J = 0 energies identical to those of the local build
  (all 4,222 levels, ZPE 10893.391847 cm⁻¹).
