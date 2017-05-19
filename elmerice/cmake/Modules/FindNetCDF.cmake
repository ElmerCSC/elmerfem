# - cmake script for finding NetCDF libraries

#  NetCDF_INCLUDE_DIR  - user modifiable choice of where netcdf headers are
#  NetCDF_LIBRARY      - user modifiable choice of where netcdf library if
#  NetCDFF_LIBRARY     - user modifiable choice of where netcdf-fortran library is
#  or NetCDF_LIBRARIES - user modifiable choice of where netcdf & netcdff libraries are

# his module returns these variables for the rest of the project to use.
#
#  NetCDF_FOUND              - True if NetCDF found including required interfaces (see below)
#  NetCDF_LIBRARIES          - All netcdf related libraries.
#  NetCDF_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If NetCDF & NetCDF_FORTRAN libraries are already defined, do nothing
IF(NetCDF_LIBRARIES AND NetCDF_INCLUDE_DIR)
   SET(NetCDF_FOUND TRUE)
   RETURN()
ELSEIF(NetCDFF_LIBRARY AND NetCDF_LIBRARY AND NetCDF_INCLUDE_DIR)
   SET(NetCDF_LIBRARIES "${NetCDFF_LIBRARY};${NetCDF_LIBRARY}" )
   SET(NetCDF_FOUND TRUE)
   RETURN()
ENDIF()

SET(NetCDF_FOUND FALSE)
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

FIND_LIBRARY(NetCDF_LIBRARY netcdf HINTS ${NETCDFLIB})
FIND_LIBRARY(NetCDFF_LIBRARY netcdff HINTS ${NETCDFFLIB})

IF (NETCDF_INCLUDE_DIR AND NetCDF_LIBRARY AND NetCDFF_LIBRARY)
  UNSET(NETCDF_FAILMSG)
  SET(NETCDFLIB_FOUND TRUE)
  SET(NetCDF_INCLUDE_DIR ${NETCDF_INCLUDE_DIR})
  SET(NetCDF_LIBRARIES "${NetCDF_LIBRARY};${NetCDFF_LIBRARY}")
ELSE()
  SET(NETCDF_FAILMSG "NetCDF libraries not found.")
ENDIF()

IF (NOT NETCDF_FAILMSG)
  SET(NetCDF_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  NETCDFINCLUDE
  NETCDFLIB
  NETCDF_FAILMSG
  NetCDF_INCLUDE_DIR
  NetCDF_LIBRARIES
  NETCDF_LIBRARY
  NETCDF_F_LIBRARY)
