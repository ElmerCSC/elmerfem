# - cmake script for finding PHDF5 libraries

#  PHDF5_INCLUDE_DIR  - user modifiable choice of where hdf5 headers are
#  PHDF5_LIBRARY      - user modifiable choice of where hdf5 library if
#  PHDF5HL_LIBRARY     - user modifiable choice of where hdf5 heigh level library is
#  or PHDF5_LIBRARIES - user modifiable choice of where hdf5 & hdf5hl libraries are

# his module returns these variables for the rest of the project to use.
#
#  PHDF5_FOUND              - True if HDF5 found including required interfaces (see below)
#  PHDF5_LIBRARIES          - All hdf5 related libraries.
#  PHDF5_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set PHDF5_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If PHDF5 & PHDF5_HL libraries are already defined, do nothing
IF(PHDF5_LIBRARIES AND PHDF5_INCLUDE_DIR)
   SET(PHDF5_FOUND TRUE)
   RETURN()
ELSEIF(PHDF5HL_LIBRARY AND PHDF5_LIBRARY AND PHDF5_INCLUDE_DIR)
   SET(PHDF5_LIBRARIES "${PHDF5HL_LIBRARY};${PHDF5_LIBRARY}" )
   SET(PHDF5_FOUND TRUE)
   RETURN()
ENDIF()

SET(PHDF5_FOUND FALSE)
SET(PHDF5INCLUDE
  "${PHDF5ROOT}/include"
  "$ENV{PHDF5ROOT}/include"
  "${PHDF5_ROOT}/include"
  "$ENV{PHDF5_ROOT}/include"
  INTERNAL
  )

FIND_PATH(PHDF5_INCLUDE_DIR
  hdf5.h 
  HINTS 
  ${PHDF5INCLUDE}
  )

SET(PHDF5LIB 
  "${PHDF5ROOT}/lib"
  "$ENV{PHDF5ROOT}/lib"
  "${PHDF5_ROOT}/lib"
  "$ENV{PHDF5_ROOT}/lib"
  INTERNAL)

SET(PHDF5HLLIB 
  "${PHDF5ROOT}/lib"
  "$ENV{PHDF5ROOT}/lib"
  "${PHDF5_ROOT}/lib"
  "$ENV{PHDF5_ROOT}/lib"
  INTERNAL)

FIND_LIBRARY(PHDF5_LIBRARY hdf5 HINTS ${PHDF5LIB})
FIND_LIBRARY(PHDF5HL_LIBRARY hdf5hl HINTS ${PHDF5HLLIB})

IF (PHDF5_INCLUDE_DIR AND PHDF5_LIBRARY AND PHDF5HL_LIBRARY)
  UNSET(PHDF5_FAILMSG)
  SET(PHDF5LIB_FOUND TRUE)
  SET(PHDF5_INCLUDE_DIR ${PHDF5_INCLUDE_DIR})
  SET(PHDF5_LIBRARIES "${PHDF5_LIBRARY};${PHDF5HL_LIBRARY}")
ELSE()
  SET(PHDF5_FAILMSG "PHDF5 libraries not found.")
ENDIF()

IF (NOT PHDF5_FAILMSG)
  SET(PHDF5_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  PHDF5INCLUDE
  PHDF5LIB
  PHDF5_FAILMSG
  PHDF5_INCLUDE_DIR
  PHDF5_LIBRARIES
  PHDF5_LIBRARY
  PHDF5HL_LIBRARY)
