# CMake toolchain file for building on taito-mic.csc.fi
#
# Author: Mikko Byckling, CSC - IT Center for Science Ltd.
# Version: 0.1

SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR x86)
SET(CMAKE_SYSTEM_VERSION 1)

# Specify the cross compilers (serial)
SET(CMAKE_C_COMPILER icc)
SET(CMAKE_Fortran_COMPILER ifort)
SET(CMAKE_CXX_COMPILER icpc)

# Specify the cross compilers (parallel)
SET(MPI_C_COMPILER mpiicc)
SET(MPI_CXX_COMPILER mpiicpc)
SET(MPI_Fortran_COMPILER mpiifort)

# Compilation flags (i.e. with optimization)
SET(CMAKE_C_FLAGS "-O3 -g" CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-O3 -g" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS "-O3 -g -fma -align all -align array64byte" CACHE STRING "")

# BLAS and LAPACK (from MKL), no threading
SET(BLAS_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_core.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_intel_sequential.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so
		   /usr/lib64/libpthread.so
		   /usr/lib64/libm.so
		   CACHE STRING "")
SET(LAPACK_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_core.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_intel_sequential.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so
		   /usr/lib64/libpthread.so
		   /usr/lib64/libm.so
		   CACHE STRING "")

SET(SCALAPACK_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_core.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_intel_sequential.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so
		   /usr/lib64/libpthread.so
		   /usr/lib64/libm.so
		   CACHE STRING "")

# Mumps
SET(MUMPSROOT /appl/opt/mumps/intel-13.1.0/intelmpi-4.1.0/4.10.0/ CACHE STRING "")
# Hypre
SET(HYPREROOT /appl/opt/hypre/intel-13.1.0/intelmpi-4.1.0/2.9.0b CACHE STRING "")

ADD_DEFINITIONS(-I$ENV{MKLROOT}/include)

