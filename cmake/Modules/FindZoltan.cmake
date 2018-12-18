# - cmake script for finding NetCDF libraries

#  ZOLTAN_INCLUDE_DIR  - user modifiable choice of where to zoltan include dir
#  ZOLTAN_LIBRARY    - user modifiable choice of where zoltan library is

# This module returns these variables for the rest of the project to use.
#
#  ZOLTAN_FOUND              - True if ZOLTAN found 
#  ZOLTAN_LIBRARY            - zoltan library is
#  ZOLTAN_INCLUDE_DIR       - zoltan include dir.

INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If ZOLTAN_LIBRARY and ZOLTAN_INCLUDE_DIR  are already defined, do nothing
IF(ZOLTAN_LIBRARY AND ZOLTAN_INCLUDE_DIR)
   SET(ZOLTAN_FOUND TRUE)
   RETURN()
ENDIF()

SET(ZOLTAN_FOUND FALSE)
SET(ZOLTANINCLUDE
  "${ZOLTANROOT}/include"
  "$ENV{ZOLTANROOT}/include"
  "${ZOLTAN_ROOT}/include"
  "$ENV{ZOLTAN_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/zoltan/include"
  INTERNAL
  )

FIND_PATH(ZOLTAN_INCLUDE_DIR
  zoltan.h
  HINTS 
  ${ZOLTANINCLUDE}
  )

SET(ZOLTANLIB 
  "${ZOLTANROOT}/lib"
  "$ENV{ZOLTANROOT}/lib"
  "${ZOLTAN_ROOT}/lib"
  "$ENV{ZOLTAN_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/zoltan/lib"
  INTERNAL)

FIND_LIBRARY(ZOLTAN_LIBRARY zoltan HINTS ${ZOLTANLIB})

IF (ZOLTAN_INCLUDE_DIR AND ZOLTAN_LIBRARY)
  UNSET(ZOLTAN_FAILMSG)
  SET(ZOLTANLIB_FOUND TRUE)
ELSE()
  SET(ZOLTAN_FAILMSG "ZOLTAN libraries not found.")
ENDIF()

IF (NOT ZOLTAN_FAILMSG)
  SET(ZOLTAN_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  ZOLTANINCLUDE
  ZOLTANLIB
  ZOLTAN_FAILMSG
  ZOLTAN_INCLUDE_DIR
  ZOLTAN_LIBRARY)
