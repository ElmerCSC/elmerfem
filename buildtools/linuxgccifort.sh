#!/bin/sh 

export CC=gcc
export CXX=g++
export F77=ifort
export FC=ifort
export ELMER_MODULES="umfpack matc post mathlibs eio hutiter fem"

sh compile.sh
