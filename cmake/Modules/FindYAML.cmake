# - cmake script for finding YAML libraries

#  YAML_INCLUDE_DIR  - user modifiable choice of where yaml headers are
#  YAML_LIBRARY      - user modifiable choice of where yaml library if
#  or YAML_LIBRARIES - user modifiable choice of where yaml libraries is 

# his module returns these variables for the rest of the project to use.
#
#  YAML_FOUND              - True if YAML found including required interfaces (see below)
#  YAML_LIBRARIES          - All yaml related libraries.
#  YAML_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set YAML_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If YAML is already defined, do nothing
IF(YAML_LIBRARIES AND YAML_INCLUDE_DIR)
   SET(YAML_FOUND TRUE)
   RETURN()
ENDIF()

SET(YAML_FOUND FALSE)
SET(YAMLINCLUDE
  "${YAMLROOT}/include"
  "$ENV{YAMLROOT}/include"
  "${YAML_ROOT}/include"
  "$ENV{YAML_ROOT}/include"
  INTERNAL
  )

FIND_PATH(YAML_INCLUDE_DIR
  libfyaml.h 
  HINTS 
  ${YAMLINCLUDE}
  )

SET(YAMLLIB 
  "${YAMLROOT}/lib"
  "$ENV{YAMLROOT}/lib"
  "${YAML_ROOT}/lib"
  "$ENV{YAML_ROOT}/lib"
  INTERNAL)

FIND_LIBRARY(YAML_LIBRARY fyaml HINTS ${YAMLLIB})

IF (YAML_INCLUDE_DIR AND YAML_LIBRARY)
  UNSET(YAML_FAILMSG)
  SET(YAMLLIB_FOUND TRUE)
  SET(YAML_INCLUDE_DIR ${YAML_INCLUDE_DIR})
  SET(YAML_LIBRARIES "${YAML_LIBRARY}")
ELSE()
  SET(YAML_FAILMSG "YAML libraries not found.")
ENDIF()

IF (NOT YAML_FAILMSG)
  SET(YAML_FOUND TRUE)
ENDIF()



MARK_AS_ADVANCED(
  YAMLINCLUDE
  YAMLLIB
  YAML_FAILMSG
  YAML_INCLUDE_DIR
  YAML_LIBRARIES
  YAML_LIBRARY
  YAMLINCLUDE
)
