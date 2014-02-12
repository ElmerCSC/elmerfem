#!/bin/bash

export CC="cc"
export CXX="CC"
export F77="f77"
export F90="f90"
export CONFFLAGS="--with-mpi-lib-dir=/opt/SUNWhpc/lib/sparcv9 --with-mpi-inc-dir=/opt/SUNWhpc/include/v9" 
export NPROCS=42

bash compile.sh
