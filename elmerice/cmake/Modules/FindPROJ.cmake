# - cmake script for finding PROJ library

#  PROJ_LIBRARY      - user modifiable choice of where proj libraries are

# his module returns these variables for the rest of the project to use.
#
#  PROJ_FOUND            - True if proj found
#  PROJ_LIBRARY          - PROJ related library

# # handle the QUIETLY and REQUIRED arguments and set PROJ_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If libraries is already defined, do nothing 
IF(PROJ_LIBRARIES)
  SET(PROJ_FOUND TRUE)
  RETURN()
ENDIF()

SET(PROJ_FOUND FALSE)

SET(PROJLIB 
  "${PROJROOT}/lib"
  "$ENV{PROJROOT}/lib"
  "${PROJ_ROOT}/lib"
  "$ENV{PROJ_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/proj/lib"
  "${PROJROOT}/lib64"
  "$ENV{PROJROOT}/lib64"
  "${PROJ_ROOT}/lib64"
  "$ENV{PROJ_ROOT}/lib64"
  "${CMAKE_SOURCE_DIR}/proj/lib64"
  INTERNAL)

FIND_LIBRARY(PROJ_LIB proj HINTS ${PROJLIB})

IF (PROJ_LIB)
  UNSET(PROJ_FAILMSG)
  SET(PROJLIB_FOUND TRUE)
  SET(PROJ_LIBRARIES ${PROJ_LIB})
ELSE()
  SET(PROJ_FAILMSG "PROJ library not found.")
ENDIF()

IF (NOT PROJ_FAILMSG)
  SET(PROJ_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(  
  PROJLIB
  PROJ_FAILMSG
  PROJ_LIBRARIES )

