# - cmake script for finding NetCDF libraries

#  MMG_INCLUDE_DIR  - user modifiable choice of where to mmg include dir
#  MMG_LIBRARY    - user modifiable choice of where mmg library is

# his module returns these variables for the rest of the project to use.
#
#  MMG_FOUND              - True if MMG found 
#  MMG_LIBRARY            - mmg library is
#  MMG_INCLUDE_DIR       - mmg include dir.

INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If MMG_LIBRARY and MMG_INCLUDE_DIR  are already defined, do nothing
IF(MMG_LIBRARY AND MMG_INCLUDE_DIR)
   SET(MMG_FOUND TRUE)
   RETURN()
ENDIF()

SET(MMG_FOUND FALSE)
SET(MMGINCLUDE
  "${MMGROOT}/include"
  "$ENV{MMGROOT}/include"
  "${MMG_ROOT}/include"
  "$ENV{MMG_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/mmg/include"
  INTERNAL
  )

FIND_PATH(MMG_INCLUDE_DIR
  mmg/libmmgf.h
  HINTS 
  ${MMGINCLUDE}
  )

SET(MMGLIB 
  "${MMGROOT}/lib"
  "$ENV{MMGROOT}/lib"
  "${MMG_ROOT}/lib"
  "$ENV{MMG_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/mmg/lib"
  INTERNAL)

FIND_LIBRARY(MMG_LIBRARY mmg HINTS ${MMGLIB})

IF (MMG_INCLUDE_DIR AND MMG_LIBRARY)
  UNSET(MMG_FAILMSG)
  SET(MMGLIB_FOUND TRUE)
ELSE()
  SET(MMG_FAILMSG "MMG libraries not found.")
ENDIF()

IF (NOT MMG_FAILMSG)
  SET(MMG_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  MMGINCLUDE
  MMGLIB
  MMG_FAILMSG
  MMG_INCLUDE_DIR
  MMG_LIBRARY)
