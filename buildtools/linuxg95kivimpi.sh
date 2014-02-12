#!/bin/sh 

export CC=gcc
export CXX=g++
export F77=g77
export FC=g95
export ELMER_MODULES="umfpack meshgen2d elmergrid matc post mathlibs eio hutiter fem"
export USE_OWN_MATHLIBS="yes"
export CONFFLAGS="--with-mpi-dir=/opt/mpich2/gnu64/"

sh compile.sh
