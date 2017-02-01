# cmake script for finding SuperLU sparse direct solver
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(NOT SuperLU_FIND_QUIETLY)
  MESSAGE(STATUS "Finding superlu from ${SuperLU_PATH}")
ENDIF()

# If SuperLU libraries are already defined, do nothing
IF(SuperLU_LIBRARIES AND SuperLU_INCLUDE_DIR)
  SET(SuperLU_FOUND TRUE)
  RETURN()
ENDIF()

SET(SuperLU_FOUND FALSE)
MESSAGE(STATUS "Finding SuperLU")

# Try to find SuperLU
FIND_PATH(SuperLU_INCLUDE_DIR 
  supermatrix.h
  PATHS 
  "${SuperLU_PATH}/include"
  )

FIND_LIBRARY(SuperLU_LIBRARIES
  NAMES
  superlu
  libsuperlu.a
  libsuperlu.so
  PATHS
  "${SuperLU_PATH}/lib"
  "${SuperLU_PATH}"
  NO_DEFAULT_PATH
  )

IF (SuperLU_INCLUDE_DIR AND SuperLU_LIBRARIES)
  UNSET(SuperLU_FAILMSG)
  SET(SuperLU_FOUND TRUE)
ELSE()
  SET(SuperLU_FAILMSG "SuperLU library not found.")
ENDIF()

IF (SuperLU_FOUND)
  IF (NOT SuperLU_FIND_QUIETLY)
    MESSAGE(STATUS "SuperLU library found.")
    MESSAGE(STATUS "SuperLU include dir: ${SuperLU_INCLUDE_DIR}")
    MESSAGE(STATUS "SuperLU libraries: ${SuperLU_LIBRARIES}")
  ENDIF()
ELSE()
  IF (SuperLU_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${SuperLU_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  SuperLU_LIBRARIES 
  SuperLU_INCLUDE_DIR 
  SuperLU_FOUND
  SuperLU_FAILMSG)
