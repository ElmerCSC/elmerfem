#!/bin/sh 

export CC=pathcc
export CXX=pathCC
export F77=pathf90
export FC=pathf90
export ELMER_MODULES="eio elmergrid meshgen2d matc mathlibs umfpack hutiter fem post"

sh compile.sh
