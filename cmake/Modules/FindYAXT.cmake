# - cmake script for finding YAXT libraries

#  YAXT_INCLUDE_DIR  - user modifiable choice of where yac headers are
#  YAXT_LIBRARY      - user modifiable choice of where yaml library if
#  or YAXT_LIBRARIES - user modifiable choice of where yaml libraries is 

# his module returns these variables for the rest of the project to use.
#
#  YAXT_FOUND              - True if YAXT found including required interfaces (see below)
#  YAXT_LIBRARIES          - All yaml related libraries.
#  YAXT_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set YAXT_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If YAXT is already defined, do nothing
IF(YAXT_LIBRARIES AND YAXT_INCLUDE_DIR)
   SET(YAXT_FOUND TRUE)
   RETURN()
ENDIF()

SET(YAXT_FOUND FALSE)
SET(YAXTINCLUDE
  "${YAXTROOT}/include"
  "$ENV{YAXTROOT}/include"
  "${YAXT_ROOT}/include"
  "$ENV{YAXT_ROOT}/include"
  INTERNAL
  )

FIND_PATH(YAXT_INCLUDE_DIR
  yaxt.h 
  HINTS 
  ${YAXTINCLUDE}
  )

SET(YAXTLIB 
  "${YAXTROOT}/lib"
  "$ENV{YAXTROOT}/lib"
  "${YAXT_ROOT}/lib"
  "$ENV{YAXT_ROOT}/lib"
  INTERNAL)

FIND_LIBRARY(YAXT_LIBRARY yaxt HINTS ${YAXTLIB})

IF (YAXT_INCLUDE_DIR AND YAXT_LIBRARY)
  UNSET(YAXT_FAILMSG)
  SET(YAXTLIB_FOUND TRUE)
  SET(YAXT_INCLUDE_DIR ${YAXT_INCLUDE_DIR})
  SET(YAXT_LIBRARIES "${YAXT_LIBRARY}")
ELSE()
  SET(YAXT_FAILMSG "YAXT libraries not found.")
ENDIF()

IF (NOT YAXT_FAILMSG)
  SET(YAXT_FOUND TRUE)
ENDIF()



MARK_AS_ADVANCED(
  YAXTINCLUDE
  YAXTLIB
  YAXT_FAILMSG
  YAXT_INCLUDE_DIR
  YAXT_LIBRARIES
  YAXT_LIBRARY
  YAXTINCLUDE
)
