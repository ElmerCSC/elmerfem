#!/bin/sh 

export CC=gcc-3.3
export CXX=g++-3.3
export F77=gfortran
export FC=gfortran
export ELMER_MODULES="umfpack meshgen2d elmergrid matc post mathlibs eio front hutiter fem post front"
export USE_OWN_MATHLIBS="yes"

sh compile.sh

