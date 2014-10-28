MESSAGE(STATUS "Finding BLACS")

SET(BLACS_FOUND FALSE)
FIND_LIBRARY(BLACS_LIBRARIES
  NAMES 
  "blacs" "blacs-pvm" "blacs-mpi" "blacs-mpich" 
  "blacs-mpich2" "blacs-openmpi" "blacs-lam"
  HINTS
  "${SCALAPACKROOT}"
  "${SCALAPACKROOT}/lib" 
  "$ENV{SCALAPACKROOT}"
  "$ENV{SCALAPACKROOT}/lib"
  "$ENV{SCALAPACK_ROOT}"
  "$ENV{SCALAPACK_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/scalapack"
  "${CMAKE_SOURCE_DIR}/scalapack/lib"
  "${BLACSROOT}/lib" 
  "$ENV{BLACSROOT}/lib"
  "$ENV{BLACS_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/scalapack/BLACS"
  "${CMAKE_SOURCE_DIR}/scalapack/BLACS/lib")

IF(BLACS_LIBRARIES)
  SET(BLACS_FOUND TRUE)
  IF (NOT BLACS_FIND_QUIETLY)
    MESSAGE(STATUS "A library with BLACS API found.")
    MESSAGE(STATUS "BLACS libraries: ${BLACS_LIBRARIES}")
  ENDIF()
ELSE()
  IF (BLACS_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "BLACS not found, needed by found SCALAPACK")
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  BLACS_FOUND
  BLACS_LIBRARIES
  )