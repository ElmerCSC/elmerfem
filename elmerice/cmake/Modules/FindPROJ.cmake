# - cmake script for finding PROJ library

#  PROJ_INCLUDE_DIR  - user modifiable choice of where proj headers are
#  PROJ_LIBRARY      - user modifiable choice of where proj libraries are

# his module returns these variables for the rest of the project to use.
#
#  PROJ_FOUND            - True if proj found
#  PROJ_LIBRARY          - PROJ related library
#  PROJ_INCLUDE_DIR      - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set PROJ_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If libraries are already defined, do nothing 
IF(PROJ_LIBRARY AND PROJ_INCLUDE_DIR)
  SET(PROJ_FOUND TRUE)
  RETURN()
ENDIF()

SET(PROJ_FOUND FALSE)
SET(PROJINCLUDE
  "${PROJROOT}/include"
  "$ENV{PROJROOT}/include"
  "${PROJ_ROOT}/include"
  "$ENV{PROJ_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/proj/include"
  INTERNAL
  )

# Try to find PROJ
FIND_PATH(PROJ_INCLUDE_DIR
  proj.mod
  HINTS 
  ${PROJINCLUDE}
  )

SET(PROJLIB 
  "${PROJROOT}/lib"
  "$ENV{PROJROOT}/lib"
  "${PROJ_ROOT}/lib"
  "$ENV{PROJ_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/proj/lib"
  INTERNAL)

FIND_LIBRARY(PROJ_LIB proj HINTS ${PROJLIB})

IF (PROJ_INCLUDE_DIR AND PROJ_LIB)
  UNSET(PROJ_FAILMSG)
  SET(PROJLIB_FOUND TRUE)
  SET(PROJ_LIBRARY ${PROJ_LIB})
ELSE()
  SET(PROJ_FAILMSG "PROJ library not found.")
ENDIF()

IF (NOT PROJ_FAILMSG)
  SET(PROJ_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(  
  PROJINCLUDE
  PROJLIB
  PROJ_FAILMSG
  PROJ_INCLUDE_DIR 
  PROJ_LIBRARY )

