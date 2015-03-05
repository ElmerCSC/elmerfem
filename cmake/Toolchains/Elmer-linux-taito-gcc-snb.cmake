# CMake toolchain file for building on taito.csc.fi
# with the GNU toolchain and Intel MKL
#
# Author: Mikko Byckling, CSC - IT Center for Science Ltd.
# Version: 0.2

# No need to cross-compile for SNB nodes on taito.csc.fi,
# i.e., do NOT set CMAKE_SYSTEM_NAME etc. here 

# Specify the cross compilers (serial)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_Fortran_COMPILER gfortran)
SET(CMAKE_CXX_COMPILER g++)

# Specify the cross compilers (parallel)
SET(MPI_C_COMPILER mpicc)
SET(MPI_CXX_COMPILER mpicxx)
SET(MPI_Fortran_COMPILER mpif90)

# Compilation flags (i.e. with optimization)
SET(CMAKE_C_FLAGS "-O2 -g" CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-O2 -g" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS "-O2 -g" CACHE STRING "")

# BLAS and LAPACK (from MKL), no threading as MPI_INIT is used
SET(BLAS_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_gf_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_core.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_sequential.so
		   /usr/lib64/libpthread.so
		   /usr/lib64/libm.so
		   CACHE STRING "")
SET(LAPACK_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_gf_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_core.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_sequential.so
		   /usr/lib64/libpthread.so
		   /usr/lib64/libm.so
		   CACHE STRING "")
SET(SCALAPACK_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_gf_lp64.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_core.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_sequential.so
		   $ENV{MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so
		   /usr/lib64/libpthread.so
		   /usr/lib64/libm.so
		   CACHE STRING "")

ADD_DEFINITIONS(-m64 -I$ENV{MKLROOT}/include -Wl,--no-as-needed)