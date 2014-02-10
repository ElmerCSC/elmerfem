#!/bin/sh 

export CC=gcc-3.3
export CXX=g++-3.3
export F77=ifort
export FC=ifort
export ELMER_MODULES="elmergrid meshgen2d umfpack matc post mathlibs eio hutiter fem front"

# get intel compiler variables in order
source /opt/intel/fc/9.0/bin/ifortvars.sh
source /opt/intel/cc/9.0/bin/iccvars.sh

sh compile.sh
