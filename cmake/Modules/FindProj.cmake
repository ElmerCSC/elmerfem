# - cmake script for finding Proj libraries

#  Proj_INCLUDE_DIR  - user modifiable choice of where yaml headers are
#  Proj_LIBRARY      - user modifiable choice of where yaml library if
#  or Proj_LIBRARIES - user modifiable choice of where yaml libraries is 

# his module returns these variables for the rest of the Project to use.
#
#  Proj_FOUND              - True if Proj found including required interfaces (see below)
#  Proj_LIBRARIES          - All yaml related libraries.
#  Proj_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set Proj_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If Proj is already defined, do nothing
IF(Proj_LIBRARIES AND Proj_INCLUDE_DIR)
   SET(Proj_FOUND TRUE)
   RETURN()
ENDIF()

SET(Proj_FOUND FALSE)
SET(ProjINCLUDE
  "${ProjROOT}/include"
  "$ENV{ProjROOT}/include"
  "${Proj_ROOT}/include"
  "$ENV{Proj_ROOT}/include"
  INTERNAL
  )

FIND_PATH(Proj_INCLUDE_DIR
  proj.h 
  HINTS 
  ${ProjINCLUDE}
  )

SET(ProjLIB 
  "${ProjROOT}/lib"
  "$ENV{ProjROOT}/lib"
  "${Proj_ROOT}/lib"
  "$ENV{Proj_ROOT}/lib"
  INTERNAL)

FIND_LIBRARY(Proj_LIBRARY proj HINTS ${ProjLIB})

IF (Proj_INCLUDE_DIR AND Proj_LIBRARY)
  UNSET(Proj_FAILMSG)
  SET(ProjLIB_FOUND TRUE)
  SET(Proj_INCLUDE_DIR ${Proj_INCLUDE_DIR})
  SET(Proj_LIBRARIES "${Proj_LIBRARY}")
ELSE()
  SET(Proj_FAILMSG "Proj libraries not found.")
ENDIF()

IF (NOT Proj_FAILMSG)
  SET(Proj_FOUND TRUE)
ENDIF()



MARK_AS_ADVANCED(
  ProjINCLUDE
  ProjLIB
  Proj_FAILMSG
  Proj_INCLUDE_DIR
  Proj_LIBRARIES
  Proj_LIBRARY
  ProjINCLUDE
)
