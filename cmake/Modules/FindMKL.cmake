# CMake script for finding MKL for Elmer (dynamic linking, Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (WITH_MPI)
  SET(SCALAPACK_NEEDED TRUE)
ELSE()
  SET(SCALAPACK_NEEDED FALSE)
ENDIF()
IF (WITH_OpenMP)
  SET(THREADS_NEEDED TRUE)
ELSE()
  SET(THREADS_NEEDED FALSE)
ENDIF()


# If MKL_LIBRARIES libraries are already defined, do nothing
IF(MKL_BLAS_LIBRARIES AND MKL_LAPACK_LIBRARIES
    AND MKL_Fortran_FLAGS)
  
  IF(NOT SCALAPACK_NEEDED OR 
         (SCALAPACK_NEEDED AND MKL_SCALAPACK_LIBRARIES))
    SET(MKL_FOUND TRUE)
    RETURN()
  ENDIF()
ENDIF()
SET(MKL_FOUND FALSE)

IF (NOT UNIX)
  MESSAGE(STATUS "Manually specify:" "MKL_BLAS_LIBRARIES " 
                                     "MKL_LAPACK_LIBRARIES "
				     "MKL_Fortran_FLAGS")
  IF(SCALAPACK_NEEDED)
    MESSAGE(STATUS "Manually specify:" "MKL_SCALAPACK_LIBRARIES")
  ENDIF()
  MESSAGE(FATAL_ERROR "Finding MKL libraries for Elmer not yet implemented for other platforms than Linux.")
ENDIF()

SET(MKLINCLUDE
  "$ENV{MKLROOT}/include"
  "${MKLROOT}/include"
  "$ENV{MKL_ROOT}/include"
  "${MKL_ROOT}/include"
  INTERNAL
  )

IF (CMAKE_SIZEOF_VOID_P MATCHES 8)
  SET(MKL_SUFFIX "_lp64")
  # 64-bit libraries
  SET(MKLLIB
    "$ENV{MKLROOT}/lib/intel64"
    "${MKLROOT}/lib/intel64"
    "$ENV{MKL_ROOT}/lib/intel64"
    "${MKL_ROOT}/lib/intel64"
    INTERNAL
    )
ELSE()
  SET(MKL_SUFFIX "")
  # 32-bit libraries
  SET(MKLLIB
    "$ENV{MKLROOT}/lib/ia32"
    "${MKLROOT}/lib/ia32"
    "$ENV{MKL_ROOT}/lib/ia32"
    "${MKL_ROOT}/lib/ia32"
    INTERNAL
    )
ENDIF()

# Set library names and compiler arguments
IF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  SET(MKL_BASENAME "gf")
  SET(MKL_THR_BASENAME "gnu")
  SET(MKL_Fortran_FLAGS "-Wl,--no-as-needed" CACHE STRING "MKL Fortran flags")
ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
 # Core libraries
 SET(MKL_BASENAME "intel")
 SET(MKL_THR_BASENAME "intel")
 SET(MKL_Fortran_FLAGS "" CACHE STRING "MKL Fortran flags")
ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
   # Core libraries
 SET(MKL_BASENAME "intel")
 SET(MKL_THR_BASENAME "pgi")
 SET(MKL_Fortran_FLAGS "" CACHE STRING "MKL Fortran flags")
ELSE()
  MESSAGE(FATAL_ERROR "Finding MKL libraries not implemented for 
                ${CMAKE_Fortran_COMPILER_ID}")
ENDIF()

# From MKL link line advisor
#GNU, seq  -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
#GNU, mt  -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm
#Intel, seq -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
#Intel, mt -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm

# Core library
SET(MKL_CORE_LIB_NAME "mkl_core" INTERNAL)
SET(MKL_NUM_LIB_NAME "mkl_${MKL_BASENAME}${MKL_SUFFIX}" INTERNAL)
# Threading library
IF (THREADS_NEEDED)
  SET(MKL_THREAD_LIB_NAME "mkl_${MKL_THR_BASENAME}_thread" INTERNAL)
ELSE()
  SET(MKL_THREAD_LIB_NAME "mkl_sequential" INTERNAL)
ENDIF()
SET(MKL_PTHREAD_LIB_NAME "pthread" INTERNAL)
SET(MKL_M_LIB_NAME "m" INTERNAL)

# Find MKL include path
FIND_PATH(MKL_INCLUDE_DIR 
  mkl.h 
  HINTS 
  ${MKLINCLUDE}
  )

# Find MKL core, BLAS and LAPACK libraries
FIND_LIBRARY(MKL_CORE_LIB ${MKL_CORE_LIB_NAME} HINTS ${MKLLIB})
FIND_LIBRARY(MKL_NUM_LIB ${MKL_NUM_LIB_NAME} HINTS ${MKLLIB})
FIND_LIBRARY(MKL_THREAD_LIB ${MKL_THREAD_LIB_NAME} HINTS ${MKLLIB})
FIND_LIBRARY(MKL_PTHREAD_LIB ${MKL_PTHREAD_LIB_NAME} HINTS ${MKLLIB})
FIND_LIBRARY(MKL_M_LIB ${MKL_M_LIB_NAME} HINTS ${MKLLIB})

IF (MKL_INCLUDE_DIR AND MKL_CORE_LIB AND 
    MKL_NUM_LIB AND MKL_THREAD_LIB AND 
    MKL_PTHREAD_LIB AND MKL_M_LIB)
  UNSET(MKL_FAILMSG)
  # SET(MKL_FOUND TRUE)
  SET(MKL_BLAS_LIBRARIES ${MKL_NUM_LIB} 
                         ${MKL_CORE_LIB} 
                         ${MKL_THREAD_LIB}
                         ${MKL_PTHREAD_LIB}
                         ${MKL_M_LIB}  CACHE FILEPATH "MKL BLAS")
  SET(MKL_LAPACK_LIBRARIES ${MKL_NUM_LIB} 
                           ${MKL_CORE_LIB} 
                           ${MKL_THREAD_LIB}
                           ${MKL_PTHREAD_LIB}
                           ${MKL_M_LIB}  CACHE FILEPATH "MKL LAPACK")
ELSE()
  SET(MKL_FAILMSG "MKL core libraries not found.")
ENDIF()

# Find BLACS and Scalapack
SET(MKL_CPARDISO_FOUND FALSE)
IF (SCALAPACK_NEEDED AND NOT MKL_FAILMSG)
  # From MKL link line advisor
  # GNU, seq:  -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm
  # GNU, mt:  -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lmkl_blacs_intelmpi_lp64 -ldl -lpthread -lm
  # Intel, seq:  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm
  # Intel, mt:  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm
  SET(MKL_BLACS_LIB_NAME "mkl_blacs_intelmpi${MKL_SUFFIX}")
  SET(MKL_SCALAPACK_LIB_NAME "mkl_scalapack${MKL_SUFFIX}")

  FIND_LIBRARY(MKL_BLACS_LIB ${MKL_BLACS_LIB_NAME} HINTS ${MKLLIB})
  FIND_LIBRARY(MKL_SCALAPACK_LIB ${MKL_SCALAPACK_LIB_NAME} HINTS ${MKLLIB})

  IF (MKL_BLACS_LIB AND MKL_SCALAPACK_LIB)
    SET(MKL_SCALAPACK_LIBRARIES ${MKL_SCALAPACK_LIB}
                                ${MKL_CORE_LIB} 
                                ${MKL_NUM_LIB} 
                                ${MKL_THREAD_LIB}
                                ${MKL_PTHREAD_LIB}
				${MKL_BLACS_LIB}
                                ${MKL_M_LIB} CACHE FILEPATH "MKL SCALAPACK libraries")
    # Attempt to find Cluster PARDISO for OpenMP compilations
    IF(THREADS_NEEDED)
      EXECUTE_PROCESS(COMMAND ${CMAKE_NM} ${MKL_NUM_LIB}
	OUTPUT_VARIABLE MKL_CPARDISO_OUTPUT
        ERROR_VARIABLE MKL_CPARDISO_ERROR)
      STRING(FIND "${MKL_CPARDISO_OUTPUT}" "cluster_sparse_solver" 
        MKL_CPARDISO_STR)
      IF("${MKL_CPARDISO_STR}" GREATER "-1" AND
	  "${MKL_CPARDISO_ERROR}" STREQUAL "")
	SET(MKL_CPARDISO_FOUND TRUE)
      ENDIF()
    ENDIF()
  ELSE()
    SET(MKL_FAILMSG "MKL BLACS and SCALAPACK libraries not found.")
  ENDIF()
ENDIF()

IF (NOT MKL_FAILMSG)
  SET(MKL_FOUND TRUE)
ENDIF()

IF (MKL_FOUND)
  IF (NOT MKL_FIND_QUIETLY)
    MESSAGE(STATUS "A library with MKL API found.")
    MESSAGE(STATUS "MKL include dir: ${MKL_INCLUDE_DIR}")
    MESSAGE(STATUS "MKL Fortran flags: ${MKL_Fortran_FLAGS}")
    MESSAGE(STATUS "MKL BLAS libraries: ${MKL_BLAS_LIBRARIES}")
    MESSAGE(STATUS "MKL LAPACK libraries: ${MKL_LAPACK_LIBRARIES}")
    IF(SCALAPACK_NEEDED)
      MESSAGE(STATUS "MKL SCALAPACK libraries: ${MKL_SCALAPACK_LIBRARIES}")
    ENDIF()
    IF(MKL_CPARDISO_FOUND)
      MESSAGE(STATUS "MKL Cluster PARDISO found")
    ENDIF()
  ENDIF()
ELSE()
  IF (MKL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${MKL_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  MKLINCLUDE
  MKLLIB
  MKL_FAILMSG
  MKL_INCLUDE_DIR
  MKL_Fortran_FLAGS
  MKL_NUM_LIB
  MKL_CORE_LIB
  MKL_THREAD_LIB
  MKL_PTHREAD_LIB
  MKL_M_LIB
  MKL_BLACS_LIB
  MKL_SCALAPACK_LIB
  MKL_BLAS_LIBRARIES
  MKL_LAPACK_LIBRARIES
  MKL_SCALAPACK_LIBRARIES
  MKL_CPARDISO_OUTPUT
  MKL_CPARDISO_STR
  MKL_CPARDISO_ERROR
  MKL_CPARDISO_FOUND
  )
