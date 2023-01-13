# - cmake script for finding FORTRANGIS library

#  FORTRANGIS_INCLUDE_DIR  - user modifiable choice of where fortrangis headers are
#  FORTRANGIS_LIBRARY      - user modifiable choice of where fortrangis libraries are

# his module returns these variables for the rest of the project to use.
#
#  FORTRANGIS_FOUND            - True if fortrangis found
#  FORTRANGIS_LIBRARY          - FORTRANGIS related library
#  FORTRANGIS_INCLUDE_DIR      - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set FORTRANGIS_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If libraries are already defined, do nothing 
IF(FORTRANGIS_LIBRARY AND FORTRANGIS_INCLUDE_DIR)
  SET(FORTRANGIS_FOUND TRUE)
  RETURN()
ENDIF()

SET(FORTRANGIS_FOUND FALSE)
SET(FORTRANGISINCLUDE
  "${FORTRANGISROOT}/include"
  "$ENV{FORTRANGISROOT}/include"
  "${FORTRANGIS_ROOT}/include"
  "$ENV{FORTRANGIS_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/fortrangis/include"
  INTERNAL
  )

# Try to find FORTRANGIS
FIND_PATH(FORTRANGIS_INCLUDE_DIR
  fortranc.mod 
  HINTS 
  ${FORTRANGISINCLUDE}
  )

SET(FORTRANGISLIB 
  "${FORTRANGISROOT}/lib"
  "$ENV{FORTRANGISROOT}/lib"
  "${FORTRANGIS_ROOT}/lib"
  "$ENV{FORTRANGIS_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/fortrangis/lib"
  INTERNAL)

FIND_LIBRARY(FORTRANGIS_LIB fortrangis HINTS ${FORTRANGISLIB})

IF (FORTRANGIS_INCLUDE_DIR AND FORTRANGIS_LIB)
  UNSET(FORTRANGIS_FAILMSG)
  SET(FORTRANGISLIB_FOUND TRUE)
  SET(FORTRANGIS_LIBRARY ${FORTRANGIS_LIB})
ELSE()
  SET(FORTRANGIS_FAILMSG "FORTRANGIS library not found.")
ENDIF()

IF (NOT FORTRANGIS_FAILMSG)
  SET(FORTRANGIS_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(  
  FORTRANGISINCLUDE
  FORTRANGISLIB
  FORTRANGIS_FAILMSG
  FORTRANGIS_INCLUDE_DIR 
  FORTRANGIS_LIBRARY )

