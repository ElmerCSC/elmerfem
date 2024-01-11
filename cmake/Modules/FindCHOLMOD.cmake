# CMake script for finding CHOLMOD
# Thomas Zwinger, CSC - IT Center for Science Ltd.
# 2024/01
#

#  CHOLMOD_INCLUDE_DIR  - user modifiable choice of where to CHOLMOD include dir
#  CHOLMOD_LIBRARY      - user modifiable choice of where CHOLMOD library is

# his module returns these variables for the rest of the project to use.
#
#  CHOLMOD_FOUND             - True if CHOLMOD found 
#  CHOLMOD_INCLUDE_DIR       - CHOLMOD include dir.
#  CHOLMOD_LIBRARIES         - needed cuda libraries




INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If CHOLMOD libraries are already defined, do nothing
IF(CHOLMOD_LIBRARIES AND CHOLMOD_INCLUDE_DIR)
   SET(CHOLMOD_FOUND TRUE)
   RETURN()
ENDIF()

SET(CHOLMOD_FOUND FALSE)
MESSAGE(STATUS "Finding CHOLMOD")

SET(CHOLMODINCLUDE
  "${CHOLMODROOT}/include"
  "$ENV{CHOLMODROOT}/include"
  "${CHOLMODROOT}/include/suitesparse"
  "$ENV{CHOLMODROOT}/include/suitesparse"
  "${CHOLMOD_ROOT}/include"
  "$ENV{CHOLMOD_ROOT}/include"
  "${CHOLMOD_ROOT}/include/suitesparse"
  "$ENV{CHOLMOD_ROOT}/include/suitesparse"
  "${CMAKE_SOURCE_DIR}/umfpack/include"
  INTERNAL
  )
# Try to find CHOLMOD
FIND_PATH(CHOLMOD_INCLUDE_DIR 
  umfpack.h
  HINTS 
  ${CHOLMODINCLUDE}
  )

SET(CHOLMODLIB 
  "${CHOLMODROOT}/lib"
  "$ENV{CHOLMODROOT}/lib"
  "${CHOLMODROOT}/lib64"
  "$ENV{CHOLMODROOT}/lib64"
  "${CHOLMOD_ROOT}/lib"
  "$ENV{CHOLMOD_ROOT}/lib"
  "${CHOLMOD_ROOT}/lib64"
  "$ENV{CHOLMOD_ROOT}/lib64"
  "${CMAKE_SOURCE_DIR}/umfpack/lib"
  INTERNAL)

FIND_LIBRARY(CHOLMOD_LIB umfpack HINTS ${CHOLMODLIB})

IF (CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIB)
  UNSET(CHOLMOD_FAILMSG)
  SET(CHOLMODLIBS_FOUND TRUE)
  SET(CHOLMOD_LIBRARIES ${CHOLMOD_LIB})
ELSE()
  SET(CHOLMOD_FAILMSG "CHOLMOD library not found.")
ENDIF()
   
IF (NOT CHOLMOD_FAILMSG)
  SET(CHOLMOD_FOUND TRUE)
ENDIF()

IF (CHOLMOD_FOUND)
  IF (NOT CHOLMOD_FIND_QUIETLY)
    MESSAGE(STATUS "A library with CHOLMOD API found.")
    MESSAGE(STATUS "CHOLMOD include dir: ${CHOLMOD_INCLUDE_DIR}")
    MESSAGE(STATUS "CHOLMOD libraries: ${CHOLMOD_LIBRARIES}")
  ENDIF()
ELSE()
  IF (CHOLMOD_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${CHOLMOD_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  CHOLMODINCLUDE
  CHOLMODLIB
  CHOLMOD_FAILMSG
  CHOLMOD_FOUND
  CHOLMOD_INCLUDE_DIR 
  CHOLMOD_LIBRARIES)
