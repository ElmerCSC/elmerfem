#!/bin/sh 

export CC=pathcc
export CXX=pathCC
export F77=pathf90
export FC=pathf90
export CONFFLAGS="--with-64bits=no"

sh compile.sh
