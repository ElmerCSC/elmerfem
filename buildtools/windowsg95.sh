#!/bin/sh 

export CC=gcc
export CXX=g++
export F77=g95
export FC=g95wrapper
export PATH="$PATH:/usr/local/lib/gcc-lib/i686-pc-cygwin/4.0.1/"

sh compile.sh

# do installer
