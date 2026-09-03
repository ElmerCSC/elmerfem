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

  # Plain search first, config package only as a fallback. Same reasoning as in
  # FindARPACK.cmake: FIND_PACKAGE(... CONFIG) executes the distribution's
  # config file, and the one Debian and Ubuntu ship include()s a targets file
  # they do not package, which is fatal however the result is tested
  # afterwards.
  INCLUDE(FindPackageHandleStandardArgs)
  MESSAGE(STATUS "Manual search of parpack")
  # Try to find PARPACK header
  SET(ARPACKINCLUDE
    "${PARPACK_ROOT}/include"
    "$ENV{PARPACK_ROOT}/include"
    "${PARPACKROOT}/include"
    "$ENV{PARPACKROOT}/include"
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
    "${PARPACK_ROOT}/lib"
    "$ENV{PARPACK_ROOT}/lib64"
    "${PARPACKROOT}/lib"
    "$ENV{PARPACKROOT}/lib64"
    "${ARPACK_ROOT}/lib"
    "$ENV{ARPACK_ROOT}/lib64"
    "${ARPACKROOT}/lib"
    "$ENV{ARPACKROOT}/lib64"
    INTERNAL
    )
  FIND_LIBRARY(PARPACK_LIBRARIES NAMES parpack parpackng parpack-ng HINTS ${ARPACKLIB})

  # PARPACK's header is frequently not shipped under its own name: arpack and
  # parpack are one package on some distributions and share an include
  # directory on others. Accept ARPACK's before concluding the plain search
  # failed, or a perfectly good libparpack sends us into the config package
  # looking for a header that was never going to be there.
  IF(PARPACK_LIBRARIES AND NOT PARPACK_INCLUDE_DIR AND ARPACK_INCLUDE_DIR)
    SET(PARPACK_INCLUDE_DIR "${ARPACK_INCLUDE_DIR}")
  ENDIF()

  IF(NOT PARPACK_LIBRARIES OR NOT PARPACK_INCLUDE_DIR)
    # Try to find with CMake config file of upstream parpack.
    FIND_PACKAGE(PARPACK CONFIG NAMES arpack arpackng arpack-ng parpack parpackng parpack-ng)
    # IF(TARGET parpack) rather than IF(PARPACK_FOUND): a config file that ran
    # to its end sets the latter even when it created nothing.
    IF(TARGET parpack)
      GET_TARGET_PROPERTY(PARPACK_INCLUDE_DIR parpack INTERFACE_INCLUDE_DIRECTORIES)
      # Most likely arpack and parpack are packed togeher (like in Arch linux)
      # or they share the same include directory even in splitted packages (parpack-dev depends on arpack-dev)
      # So in this point ARPACK_INCLUDE_DIR has to be defined or the information is loaded into the arpack
      # interface
      IF(NOT PARPACK_INCLUDE_DIR AND TARGET arpack)
        GET_TARGET_PROPERTY(PARPACK_INCLUDE_DIR arpack INTERFACE_INCLUDE_DIRECTORIES)
      ENDIF()
      GET_TARGET_PROPERTY(PARPACK_LIBRARIES parpack IMPORTED_LOCATION_RELEASE)
      # Check if a debug build type was used
      IF(NOT PARPACK_LIBRARIES)
        GET_TARGET_PROPERTY(PARPACK_LIBRARIES parpack IMPORTED_LOCATION_DEBUG)
      ENDIF()
    ENDIF()
  ENDIF(NOT PARPACK_LIBRARIES OR NOT PARPACK_INCLUDE_DIR)

ENDIF(NOT PARPACK_FOUND)

# This checks could be inadequate because this variables are not empty if nothing found
# Other option is to use the keyword REQUIRED, but this will increase cmake version to 3.18
# https://cmake.org/cmake/help/latest/command/find_library.html
# https://cmake.org/cmake/help/latest/command/find_path.html
IF (PARPACK_INCLUDE_DIR AND PARPACK_LIBRARIES)
  SET(PARPACK_FOUND TRUE)
  #The config script was not used, define the target manually
  # UNKNOWN rather than SHARED, for the same reason as in FindARPACK.cmake: a
  # SHARED imported target without IMPORTED_IMPLIB is not linkable on Windows.
  IF(NOT TARGET parpack)
    ADD_LIBRARY(parpack UNKNOWN IMPORTED)
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
