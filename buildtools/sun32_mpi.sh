#!/bin/bash

export CC=cc
export CXX=CC
export F77=f77
export F90=f90
export CONFFLAGS="--with-64bits=no --with-mpi-dir=/opt/SUNWhpc"
export NPROCS=42

bash compile.sh