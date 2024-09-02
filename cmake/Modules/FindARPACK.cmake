#CMake module to find arpack-ng libraries
#
#
#
#
#       ARPACK_INCLUDE_DIR
#       ARPACK_LIBRARIES
#       PARPACK_INCLUDE_DIR
#       PARPACK_LIBRARIES
#
#
#
#
# If the user defined arpack library locations
IF(ARPACK_LIBRARIES AND ARPACK_INCLUDE_DIR)
  SET(ARPACK_FOUND TRUE)
  MESSAGE(STATUS "Provided arpack include path ${ARPACK_INCLUDE_DIR}")
  MESSAGE(STATUS "Provided arpack libraries ${ARPACK_LIBRARIES}")
ELSE()
  SET(ARPACK_FOUND FALSE)
ENDIF()

# If MPI is found, the user also has to supply parpack libraries
IF(MPI_FOUND AND PARPACK_LIBRARIES AND PARPACK_INCLUDE_DIR)
  SET(PARPACK_FOUND TRUE)
  MESSAGE(STATUS "Provided parpack include path ${PARPACK_INCLUDE_DIR}")
  MESSAGE(STATUS "Provided parpack libraries ${PARPACK_LIBRARIES}")
  # Both libraries are defined, return
  IF(ARPACK_FOUND AND PARPACK_FOUND)
    RETURN()
  ENDIF()
ELSE()
  #Parpack is not needed and arpack is provided, return
  IF(NOT MPI_FOUND AND ARPACK_FOUND)
    RETURN()
  ELSE()
  # PARPACK_LIBRARIES or PARPACK_INCLUDE_DIR are not defined
    SET(PARPACK_FOUND FALSE)
  ENDIF()
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
    MESSAGE(STATUS "Manual search of arpack")

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
ENDIF(NOT ARPACK_FOUND)

IF(MPI_FOUND AND NOT PARPACK_FOUND)
  MESSAGE(STATUS "Finding parpack libraries")
  # Try to find with CMake config file of upstream parpack.
  FIND_PACKAGE(PARPACK CONFIG NAMES arpack arpackng arpack-ng parpack parpackng parpack-ng)
  IF(PARPACK_FOUND)
    GET_TARGET_PROPERTY(PARPACK_INCLUDE_DIR parpack INTERFACE_INCLUDE_DIRECTORIES)
    # Most likely arpack and parpack are packed togeher (like in Arch linux)
    # or they share the same include directory even in splitted packages
    IF(NOT PARPACK_INCLUDE_DIR)
      SET(PARPACK_INCLUDE_DIR ${ARPACK_INCLUDE_DIR})
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

ENDIF(MPI_FOUND AND NOT PARPACK_FOUND)


MARK_AS_ADVANCED( ARPACK_LIBRARIES ARPACK_INCLUDE_DIR ARPACK_FOUND
                  PARPACK_LIBRARIES PARPACK_INCLUDE_DIR PARPACK_FOUND )
