# - cmake script for finding XIOS libraries

#  XIOS_INCLUDE_DIR  - user modifiable choice of where xios headers are
#  XIOS_LIBRARY      - user modifiable choice of where xios library if
#  or XIOS_LIBRARIES - user modifiable choice of where xios libraries is 

# his module returns these variables for the rest of the project to use.
#
#  XIOS_FOUND              - True if XIOS found including required interfaces (see below)
#  XIOS_LIBRARIES          - All xios related libraries.
#  XIOS_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set XIOS_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If XIOS is already defined, do nothing
IF(XIOS_LIBRARIES AND XIOS_INCLUDE_DIR)
   SET(XIOS_FOUND TRUE)
   RETURN()
ENDIF()

SET(XIOS_FOUND FALSE)
SET(XIOSINCLUDE
  "${XIOSROOT}/include"
  "$ENV{XIOSROOT}/include"
  "${XIOS_ROOT}/include"
  "$ENV{XIOS_ROOT}/include"
  INTERNAL
  )
SET(XIOSINC
  "${XIOSROOT}/inc"
  "$ENV{XIOSROOT}/inc"
  "${XIOS_ROOT}/inc"
  "$ENV{XIOS_ROOT}/inc"
  INTERNAL
  )

FIND_PATH(XIOS_INCLUDE_DIR
  xios.h xios.hpp
  HINTS 
  ${XIOSINCLUDE}
  ${XIOSINC}
  )

SET(XIOSLIB 
  "${XIOSROOT}/lib"
  "$ENV{XIOSROOT}/lib"
  "${XIOS_ROOT}/lib"
  "$ENV{XIOS_ROOT}/lib"
  INTERNAL)

# Prefer shared (.so/.dll/.dylib) over static lib (.a/.lib) if both are present
IF(APPLE)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES .dylib .a)
ELSEIF(WIN32)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES .dll .lib)
ELSE()
  SET(CMAKE_FIND_LIBRARY_SUFFIXES .so .a)
ENDIF()

FIND_LIBRARY(XIOS_LIBRARY xios HINTS ${XIOSLIB})

IF (XIOS_INCLUDE_DIR AND XIOS_LIBRARY)
  GET_FILENAME_COMPONENT(XIOS_LIBRARY_EXT ${XIOS_LIBRARY} EXT)
  IF(XIOS_LIBRARY_EXT MATCHES "so|dll|dylib")
    UNSET(XIOS_FAILMSG)
    SET(XIOSLIB_FOUND TRUE)
    SET(XIOS_INCLUDE_DIR ${XIOS_INCLUDE_DIR})
    SET(XIOS_LIBRARIES "${XIOS_LIBRARY}")
  ELSE()
    SET(XIOS_FAILMSG "XIOS only found as a static library. Please rebuild XIOS as a shared library (.so/.dll/.dylib) if you want to use it.")
  ENDIF()
ELSE()
  SET(XIOS_FAILMSG "XIOS not found")
ENDIF()

IF (NOT XIOS_FAILMSG)
  SET(XIOS_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  XIOSINCLUDE
  XIOSLIB
  XIOS_FAILMSG
  XIOS_INCLUDE_DIR
  XIOS_LIBRARIES
  XIOS_LIBRARY)
