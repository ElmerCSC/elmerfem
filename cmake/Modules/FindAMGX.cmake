# CMake script for finding AMGX
# Thomas Zwinger, CSC - IT Center for Science Ltd.
# 2020/03
#

#  AMGX_INCLUDE_DIR  - user modifiable choice of where to AMGX include dir
#  AMGX_LIBRARY    - user modifiable choice of where AMGX library is

# his module returns these variables for the rest of the project to use.
#
#  AMGX_FOUND              - True if AMGX found 
#  AMGX_LIBRARY           -  AMGX library
#  AMGX_INCLUDE_DIR       - AMGX include dir.
#  CUDA_LIBRARIES         - needed cuda libraries
#  CUDA_LIBDIR             - needed cuda library directory

#INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If AMGX_LIBRARY and AMGX_INCLUDE_DIR  are already defined, do nothing
IF(AMGX_LIBRARY AND AMGX_INCLUDE_DIR AND CUDA_LIBRARIES)
   SET(AMGX_FOUND TRUE)
   RETURN()
ENDIF()

FIND_PACKAGE(CUDA)
message("Cuda libraries: " ${CUDA_LIBRARIES})

SET(AMGX_FOUND FALSE)
SET(AMGXINCLUDE
  "${AMGXROOT}/include"
  "$ENV{AMGXROOT}/include"
  "${AMGX_ROOT}/include"
  "$ENV{AMGX_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/amgx/include"
  INTERNAL
  )

FIND_PATH(AMGX_INCLUDE_DIR
  amgx_c.h
  HINTS 
  ${AMGXINCLUDE}
  )

SET(AMGXLIB 
  "${AMGXROOT}/lib"
  "$ENV{AMGXROOT}/lib"
  "${AMGX_ROOT}/lib"
  "$ENV{AMGX_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/amgx/lib"
  INTERNAL)

FIND_LIBRARY(AMGX_LIBRARY amgx HINTS ${AMGXLIB})

#SET(AMGX_LIBRARIES ${AMGXLIB} ${CUDA_LIBRARIES})

IF (AMGX_INCLUDE_DIR AND AMGX_LIBRARIES)
  UNSET(AMGX_FAILMSG)
  SET(AMGXLIB_FOUND TRUE)
ELSE()
  SET(AMGX_FAILMSG "AMGX libraries not found.")
ENDIF()

IF (NOT AMGX_FAILMSG)
  SET(AMGX_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  AMGX_FAILMSG
  AMGX_FOUND
  AMGX_INCLUDE_DIR
  AMGX_LIBRARY
  CUDA_LIBRARIES
  CUDA_LIBDIR)
