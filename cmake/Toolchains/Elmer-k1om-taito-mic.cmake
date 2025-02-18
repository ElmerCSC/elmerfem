# CMake toolchain file for building Elmer for Xeon Phi 
# on taito-mic.csc.fi
#
# Based on internet article by Dmitri Mishura at 
# https://software.intel.com/en-us/articles/cross-compilation-for-intel-xeon-phi-coprocessor-with-cmake
#
# Author: Mikko Byckling, CSC - IT Center for Science Ltd.
# Version: 0.1

SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR k1om)
SET(CMAKE_SYSTEM_VERSION 1)

# Specify the cross compilers (serial)
SET(CMAKE_C_COMPILER icc)
SET(CMAKE_Fortran_COMPILER ifort)
SET(CMAKE_CXX_COMPILER icpc)

# Specify the cross compilers (parallel)
SET(MPI_C_COMPILER mpiicc)
SET(MPI_CXX_COMPILER mpiicpc)
SET(MPI_Fortran_COMPILER mpiifort)
SET(_CMAKE_TOOLCHAIN_PREFIX x86_64-k1om-linux- CACHE STRING "")
MARK_AS_ADVANCED(_CMAKE_TOOLCHAIN_PREFIX)

# Compilation flags (i.e. -mmic with optimization)
SET(CMAKE_C_FLAGS "-mmic -O3 -g -openmp" CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-mmic -O3 -g -openmp"  CACHE STRING "")
SET(CMAKE_Fortran_FLAGS 
    "-mmic -O3 -g -openmp -align all -align array64byte"  CACHE STRING "")

# K1OM target environment, MKL, MPI
SET(CMAKE_FIND_ROOT_PATH /usr/linux-k1om-4.7
			 /usr/linux-k1om-4.7/linux-k1om
			 $ENV{MKLROOT}/lib/mic
			 $ENV{I_MPI_ROOT}/mic)

# Search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM BOTH)
# Search for libraries only in the target environment
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
# Search for file in both environments
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)

# BLAS and LAPACK (from MKL)
SET(BLAS_LIBRARIES $ENV{MKLROOT}/lib/mic/libmkl_scalapack_lp64.so
		   $ENV{MKLROOT}/lib/mic/libmkl_intel_lp64.so
		   $ENV{MKLROOT}/lib/mic/libmkl_core.so
		   $ENV{MKLROOT}/lib/mic/libmkl_intel_thread.so
		   $ENV{MKLROOT}/lib/mic/libmkl_blacs_intelmpi_lp64.so
		   /usr/linux-k1om-4.7/linux-k1om/usr/lib64/libpthread.so
		   /usr/linux-k1om-4.7/linux-k1om/usr/lib64/libm.so
		   CACHE STRING "")
SET(LAPACK_LIBRARIES $ENV{MKLROOT}/lib/mic/libmkl_scalapack_lp64.so
		   $ENV{MKLROOT}/lib/mic/libmkl_intel_lp64.so
		   $ENV{MKLROOT}/lib/mic/libmkl_core.so
		   $ENV{MKLROOT}/lib/mic/libmkl_intel_thread.so
		   $ENV{MKLROOT}/lib/mic/libmkl_blacs_intelmpi_lp64.so
		   /usr/linux-k1om-4.7/linux-k1om/usr/lib64/libpthread.so
		   /usr/linux-k1om-4.7/linux-k1om/usr/lib64/libm.so
		   CACHE STRING "")
# MKL includes
ADD_DEFINITIONS(-I$ENV{MKLROOT}/include)

# Set blas and lapack libraries as advanced settings
MARK_AS_ADVANCED(BLAS_LIBRARIES LAPACK_LIBRARIES)