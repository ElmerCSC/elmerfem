#!/bin/sh 

export CC=icc
export CXX=icc
export F77=ifort
export FC=ifort
export ELMER_MODULES="elmergrid meshgen2d umfpack matc post mathlibs eio hutiter fem front"

# get intel compiler variables in order
source /opt/intel_cc_80/bin/iccvars.sh 
source /opt/intel_fc_80/bin/ifortvars.sh

sh compile.sh
