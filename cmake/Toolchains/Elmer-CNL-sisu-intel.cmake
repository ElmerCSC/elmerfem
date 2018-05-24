# CMake toolchain file for building on sisu.csc.fi
#
# Author: Sami Ilvonen, CSC - IT Center for Science Ltd.
# Version: 0.1

SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR x86_64)
SET(CMAKE_SYSTEM_VERSION 1)

# Specify the cross compilers (serial)
SET(CMAKE_C_COMPILER cc)
SET(CMAKE_Fortran_COMPILER ftn)
SET(CMAKE_CXX_COMPILER CC)

# Specify the cross compilers (parallel)
SET(MPI_C_COMPILER cc)
SET(MPI_CXX_COMPILER CC)
SET(MPI_Fortran_COMPILER ftn)

# Compilation flags (i.e. with optimization)
SET(CMAKE_C_FLAGS "-g -O2 -fPIC" CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-g -O2 -fPIC" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS "-g -O2 -fPIC" CACHE STRING "")

# Reset the default build type flags, all flags should be set above
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "" CACHE STRING "")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "" CACHE STRING "")

# Compilation flags for native build using C compiler, used for ElmerGrid
SET(NATIVE_BUILD_FLAGS "-g -O2 -fPIC -target-cpu=sandybridge" CACHE STRING "")
# from https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
# -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
# BLAS and LAPACK (from MKL), no threading as MPI_INIT is used
SET(BLAS_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.so
                   $ENV{MKLROOT}/lib/intel64/libmkl_core.so
                   $ENV{MKLROOT}/lib/intel64/libmkl_sequential.so
                   /usr/lib64/libpthread.so
                   /usr/lib64/libm.so
                   /usr/lib64/libdl.so
                   CACHE STRING "")
SET(LAPACK_LIBRARIES $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.so
                   $ENV{MKLROOT}/lib/intel64/libmkl_core.so
                   $ENV{MKLROOT}/lib/intel64/libmkl_sequential.so
                   /usr/lib64/libpthread.so
                   /usr/lib64/libm.so
                   /usr/lib64/libdl.so
                   CACHE STRING "")

# Set the result of GFortran version test for cross compiling
SET( GFORTRAN_VERSIONTEST_RUN_RESULT 
     "0"
     CACHE STRING "Result from TRY_RUN" FORCE)

SET( GFORTRAN_VERSIONTEST_RUN_RESULT__TRYRUN_OUTPUT 
     "4.9"
     CACHE STRING "Output from TRY_RUN" FORCE)
SET ( MPIEXEC "aprun" CACHE STRING "")
SET ( MPIEXEC_NUMPROC_FLAG  "-n" CACHE STRING "")
