# CMake toolchain file for building with GNU compilers
#
# Author: Mikko Byckling, CSC - IT Center for Science Ltd.
# Version: 0.1

SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR x86_64)
SET(CMAKE_SYSTEM_VERSION 1)

# Specify the cross compilers (serial)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_Fortran_COMPILER gfortran)
SET(CMAKE_CXX_COMPILER g++)

# Specify the cross compilers (parallel)
SET(MPI_C_COMPILER mpicc)
SET(MPI_CXX_COMPILER mpic++)
SET(MPI_Fortran_COMPILER mpif90)

# Compilation flags (i.e. with optimization)
SET(CMAKE_C_FLAGS "-O3 -g -m64 -ftree-vectorize -funroll-loops" CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-O3 -g -m64 -ftree-vectorize -funroll-loops" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS "-O3 -g -m64 -ftree-vectorize -funroll-loops" CACHE STRING "")
