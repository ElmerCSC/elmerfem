#!/bin/sh 

export CC=gcc-3.3
export CXX=g++-3.3
export F77=ifc
export FC=ifc

source /opt/intel/compiler70/ia32/bin/ifcvars.sh

sh compile.sh
