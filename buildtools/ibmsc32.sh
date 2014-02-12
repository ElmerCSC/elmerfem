#!/bin/sh

export CC=xlc
export CXX=xlC
export F77=xlf
export FC=xlf90
export CONFFLAGS="--with-64bits=no"
export NPROCS=42

bash compile.sh
