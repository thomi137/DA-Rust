[![Rust CI Test](https://github.com/thomi137/DA-Rust/actions/workflows/rust.yml/badge.svg?branch=test)](https://github.com/thomi137/DA-Rust/actions/workflows/rust.yml)

# DA-Rust
Diploma thesis in physics ETH Zurich_
Some code in Rust for learning purposes

## Purpose
This repository comprises the source code used to simulate kicked Bose-Einstein condensates trapped in a lattice with harmonic potential wells. I specifically is used to examine behaviour of the governing Gross-Pitaevskii equation on specific Kick strengths. The non-linearity implies the validity of the KAM theory of classical chaos theory.

At the time it was thought that this gives rise to chaotic behaviour in the quantum realm.
I ported some of the code to rust for scientific purposes (i. e. for fun). I think most of the problems I dealt with back then
could be addressed in this modern language, although I still use LAPACK and FFTW which are second to none and worth the headache Fortran and C bring to the table.

## Installation Requirements
In order to be able to compile the file in this repository, you will need to have the following installed:

- A working C++ compiler with the standard library
- A working Fortran compiler if you are installing LAPACK and BLAS from scratch.
- Lapack and BLAS: Installation and usage docs founr [here](http://www.netlib.org/lapack/)
- The [Fastest Fourier Transform in the West] (https://www.fftw.org/) Library this can be installed via the fftw crate and compiled during build. But you still need a C compiler.
