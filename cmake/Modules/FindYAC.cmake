# - cmake script for finding YAC libraries

#  YAC_INCLUDE_DIR  - user modifiable choice of where yac headers are
#  YAC_LIBRARY      - user modifiable choice of where yac library if
#  or YAC_LIBRARIES - user modifiable choice of where yac libraries is

# his module returns these variables for the rest of the project to use.
#
#  YAC_FOUND              - True if YAC found including required interfaces (see below)
#  YAC_LIBRARIES          - All yac related libraries.
#  YAC_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set YAC_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If YAC is already defined, do nothing
IF(YAC_LIBRARIES AND YAC_INCLUDE_DIR)
   SET(YAC_FOUND TRUE)
   RETURN()
ENDIF()

SET(YAC_FOUND FALSE)
SET(YACINCLUDE
  "${YACROOT}/include"
  "$ENV{YACROOT}/include"
  "${YAC_ROOT}/include"
  "$ENV{YAC_ROOT}/include"
  INTERNAL
  )

FIND_PATH(YAC_INCLUDE_DIR
  yac.h
  HINTS
  ${YACINCLUDE}
  )

SET(YACLIB
  "${YACROOT}/lib"
  "$ENV{YACROOT}/lib"
  "${YAC_ROOT}/lib"
  "$ENV{YAC_ROOT}/lib"
  INTERNAL)

FIND_LIBRARY(YAC_LIBRARY yac HINTS ${YACLIB})

IF (YAC_INCLUDE_DIR AND YAC_LIBRARY)
  UNSET(YAC_FAILMSG)
  SET(YAC_INCLUDE_DIR ${YAC_INCLUDE_DIR})

  # List of YAC static libraries
  SET(YACLIB_NAMES
    libyac.a
    libyac_core.a
    libyac_mci.a
    libyac_mtime.a
  )

  # Determine directory of the found YAC library
  get_filename_component(YACLIB_DIR "${YAC_LIBRARY}" DIRECTORY)

  SET(YAC_LIBRARIES "")
  FOREACH(LIB ${YACLIB_NAMES})
    SET(YAC_LIB_PATH "${YACLIB_DIR}/${LIB}")
    IF (EXISTS "${YAC_LIB_PATH}")
      LIST(APPEND YAC_LIBRARIES "${YAC_LIB_PATH}")
    ELSE()
      # if any of the expected libraries is not found return with error
      SET(YAC_FAILMSG "YAC:           Expected library not found: ${YAC_LIB_PATH}.")
      SET(YACLIB_FOUND FALSE)
      RETURN()
    ENDIF()
  ENDFOREACH()
  SET(YACLIB_FOUND TRUE)
ELSE()
  SET(YAC_FAILMSG "YAC libraries not found.")
ENDIF()

IF (NOT YAC_FAILMSG)
  SET(YAC_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  YACINCLUDE
  YACLIB
  YAC_FAILMSG
  YAC_INCLUDE_DIR
  YAC_LIBRARIES
  YAC_LIBRARY
)
