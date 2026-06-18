# - cmake script for finding PDAF libraries
#
# this module returns these variables for the rest of the project to use.
#
#  PDAF_FOUND                - True if PDAF found including required interfaces (see below)
#  PDAF_LIBRARIES            - PDAF related libraries.
#  PDAF_INCLUDE_DIR	     - PDAF include directory
#
# # handle the QUIETLY and REQUIRED arguments and set PDAF_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

IF(PDAF_LIBRARIES AND PDAF_INCLUDE_DIR)
   SET(PDAF_FOUND TRUE)
   RETURN()
ENDIF()

SET(PDAF_FOUND FALSE)

## PDAF inc dir
SET(PDAFINC
   "${PDAFROOT}/inc"
   "$ENV{PDAFROOT}/inc"
   "${PDAF_ROOT}/inc"
   "$ENV{PDAF_ROOT}/inc"
   INTERNAL)

FIND_PATH(PDAF_INCLUDE_DIR pdaf.mod HINTS ${PDAFINC})

## PDAF lib
SET(PDAFLIB
  "${PDAFROOT}/lib"
  "$ENV{PDAFROOT}/lib"
  "${PDAF_ROOT}/lib"
  "$ENV{PDAF_ROOT}/lib"
  INTERNAL)

FIND_LIBRARY(PDAF_LIBRARY pdaf-d HINTS ${PDAFLIB})

##
IF (PDAF_INCLUDE_DIR AND PDAF_LIBRARY)
   UNSET(PDAF_FAILMSG)
   SET(PDAF_FOUND TRUE)
   SET(PDAF_INCLUDE_DIR ${PDAF_INCLUDE_DIR})
   SET(PDAF_LIBRARIES ${PDAF_LIBRARY})
ELSE()
	SET(PDAF_FAILMSG "PDAF not found.")
ENDIF()

