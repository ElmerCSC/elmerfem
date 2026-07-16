# - cmake script for finding NETCDF libraries

#  NETCDF_INCLUDE_DIR  - user modifiable choice of where netcdf headers are
#  NETCDF_LIBRARY      - user modifiable choice of where netcdf library if
#  NETCDFF_LIBRARY     - user modifiable choice of where netcdf-fortran library is
#  or NETCDF_LIBRARIES - user modifiable choice of where netcdf & netcdff libraries are

# his module returns these variables for the rest of the project to use.
#
#  NETCDF_FOUND              - True if NETCDF found including required interfaces (see below)
#  NETCDF_LIBRARIES          - All netcdf related libraries.
#  NETCDF_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If NETCDF & NETCDF_FORTRAN libraries are already defined, do nothing
IF(NETCDF_LIBRARIES AND NETCDF_INCLUDE_DIR)
   SET(NETCDF_FOUND TRUE)
   RETURN()
ELSEIF(NETCDFF_LIBRARY AND NETCDF_LIBRARY AND NETCDF_INCLUDE_DIR)
   SET(NETCDF_LIBRARIES "${NETCDFF_LIBRARY};${NETCDF_LIBRARY}" )
   SET(NETCDF_FOUND TRUE)
   RETURN()
ENDIF()

SET(NETCDF_FOUND FALSE)
SET(NETCDFINCLUDE
  "${NETCDFROOT}/include"
  "$ENV{NETCDFROOT}/include"
  "${NETCDF_ROOT}/include"
  "$ENV{NETCDF_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/netcdf/include"
  INTERNAL
  )

FIND_PATH(NETCDF_INCLUDE_DIR
  netcdf.h 
  HINTS 
  ${NETCDFINCLUDE}
  )

SET(NETCDFLIB 
  "${NETCDFROOT}/lib"
  "$ENV{NETCDFROOT}/lib"
  "${NETCDF_ROOT}/lib"
  "$ENV{NETCDF_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/netcdf/lib"
  INTERNAL)

SET(NETCDFFLIB 
  "${NETCDFROOT}/lib"
  "$ENV{NETCDFROOT}/lib"
  "${NETCDF_ROOT}/lib"
  "$ENV{NETCDF_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/netcdf/lib"
  INTERNAL)

FIND_LIBRARY(NETCDF_LIBRARY netcdf HINTS ${NETCDFLIB})
FIND_LIBRARY(NETCDFF_LIBRARY netcdff HINTS ${NETCDFFLIB})

IF (NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY AND NETCDFF_LIBRARY)
  UNSET(NETCDF_FAILMSG)
  SET(NETCDFLIB_FOUND TRUE)
  SET(NETCDF_INCLUDE_DIR ${NETCDF_INCLUDE_DIR})
  SET(NETCDF_LIBRARIES "${NETCDFF_LIBRARY};${NETCDF_LIBRARY}")
ELSE()
  SET(NETCDF_FAILMSG "NETCDF libraries not found.")
ENDIF()

IF (NOT NETCDF_FAILMSG)
  SET(NETCDF_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  NETCDFINCLUDE
  NETCDFLIB
  NETCDF_FAILMSG
  NETCDF_INCLUDE_DIR
  NETCDF_LIBRARIES
  NETCDF_LIBRARY
  NETCDF_F_LIBRARY)
