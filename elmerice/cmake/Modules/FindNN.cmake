# - cmake script for finding NN library

#  NN_INCLUDE_DIR  - user modifiable choice of where nn headers are
#  NN_LIBRARY      - user modifiable choice of where nn libraries are

# his module returns these variables for the rest of the project to use.
#
#  NN_FOUND            - True if nn found
#  NN_LIBRARY          - NN related library
#  NN_INCLUDE_DIR      - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set NN_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If libraries are already defined, do nothing 
IF(NN_LIBRARY AND NN_INCLUDE_DIR)
  SET(NN_FOUND TRUE)
  RETURN()
ENDIF()

SET(NN_FOUND FALSE)
SET(NNINCLUDE
  "${NNROOT}/"
  "$ENV{NNROOT}/"
  "${NN_ROOT}/"
  "$ENV{NN_ROOT}/"
  "${CMAKE_SOURCE_DIR}/nn-c/"
  INTERNAL
  )

# Try to find NN
FIND_PATH(NN_INCLUDE_DIR
  nan.h 
  HINTS 
  ${NNINCLUDE}
  )

SET(NNLIB 
  "${NNROOT}/lib"
  "$ENV{NNROOT}/lib"
  "${NN_ROOT}/lib"
  "$ENV{NN_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/nn-c/lib"
  INTERNAL)

FIND_LIBRARY(NN_LIB nn HINTS ${NNLIB})

IF (NN_INCLUDE_DIR AND NN_LIB)
  UNSET(NN_FAILMSG)
  SET(NNLIB_FOUND TRUE)
  SET(NN_LIBRARY ${NN_LIB})
ELSE()
  SET(NN_FAILMSG "NN library not found.")
ENDIF()

IF (NOT NN_FAILMSG)
  SET(NN_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(  
  NNINCLUDE
  NNLIB
  NN_FAILMSG
  NN_INCLUDE_DIR 
  NN_LIBRARY )

