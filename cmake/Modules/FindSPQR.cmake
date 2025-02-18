# CMake script for finding SPQR

#  SPQR_INCLUDE_DIR  - user modifiable choice of where to SPQR include dir
#  SPQR_LIBRARY      - user modifiable choice of where SPQR library is

# his module returns these variables for the rest of the project to use.
#
#  SPQR_FOUND             - True if SPQR found
#  SPQR_INCLUDE_DIR       - SPQR include dir.
#  SPQR_LIBRARIES         - needed SPQR libraries

# If SPQR libraries are already defined, do nothing
IF(SPQR_LIBRARIES AND SPQR_INCLUDE_DIR)
   SET(SPQR_FOUND TRUE)
   RETURN()
ENDIF()

# Try to find with CMake config file of upstream UMFPACK.
FIND_PACKAGE(SPQR CONFIG)

IF(SPQR_LIBRARIES AND SPQR_INCLUDE_DIR)
  RETURN()
ENDIF()

# Fall back to manual search
INCLUDE(FindPackageHandleStandardArgs)

SET(SPQR_FOUND FALSE)
MESSAGE(STATUS "Finding SPQR")

SET(SPQRINCLUDE
  "${SPQRROOT}/include"
  "$ENV{SPQRROOT}/include"
  "${SPQRROOT}/include/suitesparse"
  "$ENV{SPQRROOT}/include/suitesparse"
  "${SPQR_ROOT}/include"
  "$ENV{SPQR_ROOT}/include"
  "${SPQR_ROOT}/include/suitesparse"
  "$ENV{SPQR_ROOT}/include/suitesparse"
  "${CMAKE_SOURCE_DIR}/spqr/include"
  INTERNAL
  )
# Try to find SPQR
FIND_PATH(SPQR_INCLUDE_DIR
  SuiteSparseQR_C.h
  HINTS
  ${SPQRINCLUDE}
  )

SET(SPQRLIB
  "${SPQRROOT}/lib"
  "$ENV{SPQRROOT}/lib"
  "${SPQRROOT}/lib64"
  "$ENV{SPQRROOT}/lib64"
  "${SPQR_ROOT}/lib"
  "$ENV{SPQR_ROOT}/lib"
  "${SPQR_ROOT}/lib64"
  "$ENV{SPQR_ROOT}/lib64"
  "${CMAKE_SOURCE_DIR}/spqr/lib"
  INTERNAL)

FIND_LIBRARY(SPQR_LIB spqr HINTS ${SPQRLIB})

IF (SPQR_INCLUDE_DIR AND SPQR_LIB)
  UNSET(SPQR_FAILMSG)
  SET(SPQRLIBS_FOUND TRUE)
  SET(SPQR_LIBRARIES ${SPQR_LIB})
ELSE()
  SET(SPQR_FAILMSG "SPQR library not found.")
ENDIF()

IF (NOT SPQR_FAILMSG)
  SET(SPQR_FOUND TRUE)
ENDIF()

IF (SPQR_FOUND)
  IF (NOT SPQR_FIND_QUIETLY)
    MESSAGE(STATUS "A library with SPQR API found.")
    MESSAGE(STATUS "SPQR include dir: ${SPQR_INCLUDE_DIR}")
    MESSAGE(STATUS "SPQR libraries: ${SPQR_LIBRARIES}")
  ENDIF()
ELSE()
  IF (SPQR_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${SPQR_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  SPQRINCLUDE
  SPQRLIB
  SPQR_FAILMSG
  SPQR_FOUND
  SPQR_INCLUDE_DIR
  SPQR_LIBRARIES)
