# CMake script for finding UMFPACK
# Thomas Zwinger, CSC - IT Center for Science Ltd.
# 2024/01
#

#  UMFPACK_INCLUDE_DIR  - user modifiable choice of where to UMFPACK include dir
#  UMFPACK_LIBRARY      - user modifiable choice of where UMFPACK library is

# his module returns these variables for the rest of the project to use.
#
#  UMFPACK_FOUND             - True if UMFPACK found 
#  UMFPACK_INCLUDE_DIR       - UMFPACK include dir.
#  UMFPACK_LIBRARIES         - needed cuda libraries




INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If UMFPACK libraries are already defined, do nothing
IF(UMFPACK_LIBRARIES AND UMFPACK_INCLUDE_DIR)
   SET(UMFPACK_FOUND TRUE)
   RETURN()
ENDIF()

SET(UMFPACK_FOUND FALSE)
MESSAGE(STATUS "Finding UMFPACK")

SET(UMFPACKINCLUDE
  "${UMFPACKROOT}/include"
  "$ENV{UMFPACKROOT}/include"
  "${UMFPACKROOT}/include/suitesparse"
  "$ENV{UMFPACKROOT}/include/suitesparse"
  "${UMFPACK_ROOT}/include"
  "$ENV{UMFPACK_ROOT}/include"
  "${UMFPACK_ROOT}/include/suitesparse"
  "$ENV{UMFPACK_ROOT}/include/suitesparse"
  "${CMAKE_SOURCE_DIR}/umfpack/include"
  INTERNAL
  )
# Try to find UMFPACK
FIND_PATH(UMFPACK_INCLUDE_DIR 
  umfpack.h
  HINTS 
  ${UMFPACKINCLUDE}
  )

SET(UMFPACKLIB 
  "${UMFPACKROOT}/lib"
  "$ENV{UMFPACKROOT}/lib"
  "${UMFPACKROOT}/lib64"
  "$ENV{UMFPACKROOT}/lib64"
  "${UMFPACK_ROOT}/lib"
  "$ENV{UMFPACK_ROOT}/lib"
  "${UMFPACK_ROOT}/lib64"
  "$ENV{UMFPACK_ROOT}/lib64"
  "${CMAKE_SOURCE_DIR}/umfpack/lib"
  INTERNAL)

FIND_LIBRARY(UMFPACK_LIB umfpack HINTS ${UMFPACKLIB})

IF (UMFPACK_INCLUDE_DIR AND UMFPACK_LIB)
  UNSET(UMFPACK_FAILMSG)
  SET(UMFPACKLIBS_FOUND TRUE)
  SET(UMFPACK_LIBRARIES ${UMFPACK_LIB})
ELSE()
  SET(UMFPACK_FAILMSG "UMFPACK library not found.")
ENDIF()
   
IF (NOT UMFPACK_FAILMSG)
  SET(UMFPACK_FOUND TRUE)
ENDIF()

IF (UMFPACK_FOUND)
  IF (NOT UMFPACK_FIND_QUIETLY)
    MESSAGE(STATUS "A library with UMFPACK API found.")
    MESSAGE(STATUS "UMFPACK include dir: ${UMFPACK_INCLUDE_DIR}")
    MESSAGE(STATUS "UMFPACK libraries: ${UMFPACK_LIBRARIES}")
  ENDIF()
ELSE()
  IF (UMFPACK_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${UMFPACK_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  UMFPACKINCLUDE
  UMFPACKLIB
  UMFPACK_FAILMSG
  UMFPACK_FOUND
  UMFPACK_INCLUDE_DIR 
  UMFPACK_LIBRARIES)