#CMake module to find arpack-ng libraries
#
#
#
#
#       PARPACK_INCLUDE_DIR
#       PARPACK_LIBRARIES
#
#
#
#
# The user defined parpack locations
IF(PARPACK_LIBRARIES AND PARPACK_INCLUDE_DIR)
  SET(PARPACK_FOUND TRUE)
  MESSAGE(STATUS "Provided parpack include path ${PARPACK_INCLUDE_DIR}")
  MESSAGE(STATUS "Provided parpack libraries ${PARPACK_LIBRARIES}")
ELSE()
  # PARPACK_LIBRARIES or PARPACK_INCLUDE_DIR are not defined
  SET(PARPACK_FOUND FALSE)
ENDIF()

IF(NOT PARPACK_FOUND)
  MESSAGE(STATUS "Finding parpack libraries")
  # Try to find with CMake config file of upstream parpack.
  FIND_PACKAGE(PARPACK CONFIG NAMES arpack arpackng arpack-ng parpack parpackng parpack-ng)
  IF(PARPACK_FOUND)
    GET_TARGET_PROPERTY(PARPACK_INCLUDE_DIR parpack INTERFACE_INCLUDE_DIRECTORIES)
    # Most likely arpack and parpack are packed togeher (like in Arch linux)
    # or they share the same include directory even in splitted packages (parpack-dev depends on arpack-dev)
    # So in this point ARPACK_INCLUDE_DIR has to be defined or the information is loaded into the arpack
    # interface
    IF(NOT PARPACK_INCLUDE_DIR)
      GET_TARGET_PROPERTY(PARPACK_INCLUDE_DIR arpack INTERFACE_INCLUDE_DIRECTORIES)
    ENDIF()
    GET_TARGET_PROPERTY(PARPACK_LIBRARIES parpack IMPORTED_LOCATION_RELEASE)
    # Check if a debug build type was used
    IF(NOT PARPACK_LIBRARIES)
      GET_TARGET_PROPERTY(PARPACK_LIBRARIES parpack IMPORTED_LOCATION_DEBUG)
    ENDIF()
  ENDIF()
  # There is no parpack-config script or something went wrong with the script
  # Fall back to manual search
  IF(NOT PARPACK_LIBRARIES OR NOT PARPACK_INCLUDE_DIR)
    INCLUDE(FindPackageHandleStandardArgs)
    MESSAGE(STATUS "Manual search of parpack")
    # Try to find PARPACK header
    SET(ARPACKINCLUDE
      "${ARPACK_ROOT}/include"
      "$ENV{ARPACK_ROOT}/include"
      "${ARPACKROOT}/include"
      "$ENV{ARPACKROOT}/include"
      INTERNAL
      )
    FIND_PATH(PARPACK_INCLUDE_DIR NAMES parpack.h parpackng.h parpack-ng.h
      HINTS ${ARPACKINCLUDE} PATH_SUFFIXES arpack arpackng arpack-ng parpack parpackng parpack-ng)
    # Try to find PARPACK libraries
    SET(ARPACKLIB
      "${ARPACK_ROOT}/lib"
      "$ENV{ARPACK_ROOT}/lib64"
      "${ARPACKROOT}/lib"
      "$ENV{ARPACKROOT}/lib64"
      INTERNAL
      )
    FIND_LIBRARY(PARPACK_LIBRARIES NAMES parpack parpackng parpack-ng HINTS ${ARPACKLIB})
  ENDIF(NOT PARPACK_LIBRARIES OR NOT PARPACK_INCLUDE_DIRR)

ENDIF(NOT PARPACK_FOUND)

# This checks could be inadequate because this variables are not empty if nothing found
# Other option is to use the keyword REQUIRED, but this will increase cmake version to 3.18
# https://cmake.org/cmake/help/latest/command/find_library.html
# https://cmake.org/cmake/help/latest/command/find_path.html
IF (PARPACK_INCLUDE_DIR AND PARPACK_LIBRARIES)
  SET(PARPACK_FOUND TRUE)
  #The config script was not used, define the target manually
  IF(NOT TARGET parpack)
    ADD_LIBRARY(parpack SHARED IMPORTED)
    SET_PROPERTY(TARGET parpack PROPERTY IMPORTED_LOCATION ${PARPACK_LIBRARIES})
    SET_PROPERTY(TARGET parpack PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PARPACK_INCLUDE_DIR})
  ENDIF()
  IF (NOT PARPACK_FIND_QUIETLY)
    MESSAGE(STATUS "A library with PARPACK API found.")
    MESSAGE(STATUS "PARPACK include dir: ${PARPACK_INCLUDE_DIR}")
    MESSAGE(STATUS "PARPACK libraries: ${PARPACK_LIBRARIES}")
  ENDIF()
ELSE()
  SET(PARPACK_FOUND FALSE)
  IF (PARPACK_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "PARPACK library not found.")
  ENDIF()
ENDIF(PARPACK_INCLUDE_DIR AND PARPACK_LIBRARIES)

MARK_AS_ADVANCED( PARPACK_LIBRARIES PARPACK_INCLUDE_DIR PARPACK_FOUND )
