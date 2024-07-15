# toolchain file for gcc8
#
# Author: Thomas Zwinger , CSC - IT Center for Science Ltd.
# Version: 1.0

# No need to cross-compile 
# i.e., do NOT set CMAKE_SYSTEM_NAME etc. here
#SET(CMAKE_SYSTEM_NAME Linux)
#SET(CMAKE_SYSTEM_PROCESSOR x86_64)
#SET(CMAKE_SYSTEM_VERSION 1) 

# Specify the compilers (serial)
SET(CMAKE_C_COMPILER gcc CACHE STRING "")
SET(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")
SET(CMAKE_CXX_COMPILER g++ CACHE STRING "")

# Specify the compilers (parallel)

# if using whole intel 18 suite
SET(MPI_C_COMPILER mpicc CACHE STRING "")
SET(MPI_CXX_COMPILER mpicxx CACHE STRING "")
SET(MPI_Fortran_COMPILER mpif90 CACHE STRING "")
