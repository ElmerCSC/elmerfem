#CMake module to find arpack-ng libraries
#
#
#
#
#       ARPACK_INCLUDE_DIR
#       ARPACK_LIBRARIES
#
#
#
#
# The user defined arpack locations
IF(ARPACK_LIBRARIES AND ARPACK_INCLUDE_DIR)
  # The user defined arpack locations
  SET(ARPACK_FOUND TRUE)
  MESSAGE(STATUS "Provided arpack include path ${ARPACK_INCLUDE_DIR}")
  MESSAGE(STATUS "Provided arpack libraries ${ARPACK_LIBRARIES}")
ELSE()
  # One location is missing
  SET(ARPACK_FOUND FALSE)
ENDIF()

IF(NOT ARPACK_FOUND)
  MESSAGE(STATUS "Finding arpack libraries")
  # Try to find with CMake config file of upstream arpack.
  FIND_PACKAGE(ARPACK CONFIG NAMES arpack arpackng arpack-ng)
  IF(ARPACK_FOUND)
    GET_TARGET_PROPERTY(ARPACK_INCLUDE_DIR arpack INTERFACE_INCLUDE_DIRECTORIES)
    GET_TARGET_PROPERTY(ARPACK_LIBRARIES arpack IMPORTED_LOCATION_RELEASE)
    # Check if a debug build type was used
    IF(NOT ARPACK_LIBRARIES)
      GET_TARGET_PROPERTY(ARPACK_LIBRARIES arpack IMPORTED_LOCATION_DEBUG)
    ENDIF()
  ENDIF()

  # There is no arpack-config script or something went wrong with the script
  # Fall back to manual search
  IF(NOT ARPACK_LIBRARIES OR NOT ARPACK_INCLUDE_DIR)
    # Fall back to manual search
    INCLUDE(FindPackageHandleStandardArgs)
    MESSAGE(STATUS "Searching for arpack library")

    # Try to find ARPACK header
    SET(ARPACKINCLUDE
      "${ARPACK_ROOT}/include"
      "$ENV{ARPACK_ROOT}/include"
      "${ARPACKROOT}/include"
      "$ENV{ARPACKROOT}/include"
      INTERNAL
      )
    FIND_PATH(ARPACK_INCLUDE_DIR NAMES arpack.h arpackng.h arpack-ng.h
      HINTS ${ARPACKINCLUDE} PATH_SUFFIXES arpack arpackng arpack-ng)

    # Try to find ARPACK libraries
    SET(ARPACKLIB
      "${ARPACK_ROOT}/lib"
      "$ENV{ARPACK_ROOT}/lib64"
      "${ARPACKROOT}/lib"
      "$ENV{ARPACKROOT}/lib64"
      INTERNAL
      )
    FIND_LIBRARY(ARPACK_LIBRARIES NAMES arpack arpackng arpack-ng HINTS ${ARPACKLIB})

  ENDIF(NOT ARPACK_LIBRARIES OR NOT ARPACK_INCLUDE_DIR)

ENDIF(NOT ARPACK_FOUND)

# This checks could be inadequate because this variables are not empty if nothing found
# Other option is to use the keyword REQUIRED, but this will increase cmake version to 3.18
# https://cmake.org/cmake/help/latest/command/find_library.html
# https://cmake.org/cmake/help/latest/command/find_path.html
IF (ARPACK_INCLUDE_DIR AND ARPACK_LIBRARIES)
  SET(ARPACK_FOUND TRUE)
  #The config script was not used, define the target manually
  IF(NOT TARGET arpack)
    ADD_LIBRARY(arpack SHARED IMPORTED)
    SET_PROPERTY(TARGET arpack PROPERTY IMPORTED_LOCATION ${ARPACK_LIBRARIES})
    SET_PROPERTY(TARGET arpack PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${ARPACK_INCLUDE_DIR})
  ENDIF()
  IF (NOT ARPACK_FIND_QUIETLY)
    MESSAGE(STATUS "A library with ARPACK API found.")
    MESSAGE(STATUS "ARPACK include dir: ${ARPACK_INCLUDE_DIR}")
    MESSAGE(STATUS "ARPACK libraries: ${ARPACK_LIBRARIES}")
  ENDIF()
ELSE()
  SET(ARPACK_FOUND FALSE)
  IF (ARPACK_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "ARPACK library not found.")
  ENDIF()
ENDIF(ARPACK_INCLUDE_DIR AND ARPACK_LIBRARIES)





MARK_AS_ADVANCED( ARPACK_LIBRARIES ARPACK_INCLUDE_DIR ARPACK_FOUND )
