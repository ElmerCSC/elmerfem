# - cmake script for finding CSA library
#
#  CSA_INCLUDE_DIR  - user modifiable choice of where csa headers are
#  CSA_LIBRARY      - user modifiable choice of where csa library is

# his module returns these variables for the rest of the project to use.
#
#  CSA_FOUND            - True if csa found
#  CSA_LIBRARY          - CSA related library
#  CSA_INCLUDE_DIR      - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set CSA_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If libraries are already defined, do nothing 
IF(CSA_LIBRARY AND CSA_INCLUDE_DIR)
  SET(CSA_FOUND TRUE)
  RETURN()
ENDIF()

SET(CSA_FOUND FALSE)
SET(CSAINCLUDE
  "${CSAROOT}/"
  "$ENV{CSAROOT}/"
  "${CSA_ROOT}/"
  "$ENV{CSA_ROOT}/"
  "${CMAKE_SOURCE_DIR}/csa-c/"
  INTERNAL
  )

# Try to find CSA
FIND_PATH(CSA_INCLUDE_DIR
  csa.h 
  HINTS 
  ${CSAINCLUDE}
  )

SET(CSALIB 
  "${CSAROOT}/lib"
  "$ENV{CSAROOT}/lib"
  "${CSA_ROOT}/lib"
  "$ENV{CSA_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/csa-c/lib"
  INTERNAL)

FIND_LIBRARY(CSA_LIB csa HINTS ${CSALIB})

IF (CSA_INCLUDE_DIR AND CSA_LIB)
  UNSET(CSA_FAILMSG)
  SET(CSALIB_FOUND TRUE)
  SET(CSA_LIBRARY ${CSA_LIB})
ELSE()
  SET(CSA_FAILMSG "CSA library not found.")
ENDIF()

IF (NOT CSA_FAILMSG)
  SET(CSA_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(  
  CSAINCLUDE
  CSALIB
  CSA_FAILMSG
  CSA_INCLUDE_DIR 
  CSA_LIBRARY )

