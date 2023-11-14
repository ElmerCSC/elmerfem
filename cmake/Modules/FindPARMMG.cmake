# - cmake script for finding NetCDF libraries

#  MMG_INCLUDE_DIR  - user modifiable choice of where to mmg include dir
#  MMG_LIBRARY    - user modifiable choice of where mmg library is

# his module returns these variables for the rest of the project to use.
#
#  PARMMG_FOUND              - True if MMG found 
#  PARMMG_LIBRARY            - mmg library is
#  PARMMG_INCLUDE_DIR       - mmg include dir.

INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If MMG_LIBRARY and MMG_INCLUDE_DIR  are already defined, do nothing
IF(PARMMG_LIBRARY AND PARMMG_INCLUDE_DIR)
   SET(PARMMG_FOUND FALSE)
   RETURN()
ENDIF()

SET(PARMMG_FOUND FALSE)
SET(PARMMGINCLUDE
  "${PARMMGROOT}/include"
  "$ENV{PARMMGROOT}/include"
  "${PARMMG_ROOT}/include"
  "$ENV{PARMMG_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/parmmg/include"
  INTERNAL
  )

FIND_PATH(PARMMG_INCLUDE_DIR
  parmmg/libparmmgf.h
  HINTS 
  ${PARMMGINCLUDE}
  )

SET(PARMMGLIB 
  "${PARMMGROOT}/lib"
  "$ENV{PARMMGROOT}/lib"
  "${PARMMG_ROOT}/lib"
  "$ENV{PARMMG_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/parmmg/lib"
  INTERNAL)

FIND_LIBRARY(PARMMG_LIBRARY parmmg HINTS ${PARMMGLIB})

IF (PARMMG_INCLUDE_DIR AND PARMMG_LIBRARY)
  UNSET(PARMMG_FAILMSG)
  SET(PARMMGLIB_FOUND TRUE)
  SET(PARMMG_FOUND TRUE)
ELSE()
  SET(PARMMG_FAILMSG "ParMMG libraries not found.")
ENDIF()

IF (NOT PARMMG_FAILMSG)
  SET(PARMMG_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  PARMMGINCLUDE
  PARMMGLIB
  PARMMG_FAILMSG
  PARMMG_INCLUDE_DIR
  PARMMG_LIBRARY)
