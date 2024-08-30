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
  # ARPACK_LIBRARIES or ARPACK_INCLUDE_DIR are not defined
    SET(PARPACK_FOUND FALSE)
  ENDIF()
ENDIF()

IF(NOT ARPACK_FOUND)
  MESSAGE(STATUS "Finding arpack libraries")
  # Try to find with CMake config file of upstream arpack.
  FIND_PACKAGE(ARPACK CONFIG NAMES arpack arpackng arpack-ng)
  IF(ARPACK_LIBRARIES AND ARPACK_INCLUDE_DIR)
    RETURN()
  ENDIF()

  # Fall back to manual search
  INCLUDE(FindPackageHandleStandardArgs)
  MESSAGE(STATUS "Finding arpack")

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
  FIND_LIBRARY(ARPACK_LIB NAMES arpack arpackng arpack-ng HINTS ${ARPACKLIB})

  # This checks may be inadequate because this variables are not empty if nothing found
# Other options is to use the keyword REQUIRED, but this will increase cmake version to 3.18
# https://cmake.org/cmake/help/latest/command/find_library.html
# https://cmake.org/cmake/help/latest/command/find_path.html
  IF (ARPACK_INCLUDE_DIR AND ARPACK_LIB)
    SET(ARPACKLIBS_FOUND TRUE)
    SET(ARPACK_FOUND TRUE)
    SET(ARPACK_LIBRARIES ${ARPACK_LIB})
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
  ENDIF()

ENDIF()

IF(MPI_FOUND AND NOT PARPACK_FOUND)
  MESSAGE(STATUS "Finding parpack libraries")
  FIND_PACKAGE(PARPACK CONFIG NAMES arpack arpackng arpack-ng parpack parpackng parpack-ng)
  IF(PARPACK_LIBRARIES AND PARPACK_INCLUDE_DIR)
    RETURN()
  ENDIF()

  # Fall back to manual search
  INCLUDE(FindPackageHandleStandardArgs)

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
  FIND_LIBRARY(PARPACK_LIB NAMES parpack parpackng parpack-ng HINTS ${ARPACKLIB})

# This checks may be inadequate because this variables are not empty if nothing found
# Other options is to use the keyword REQUIRED, but this will increase cmake version to 3.18
# https://cmake.org/cmake/help/latest/command/find_library.html
# https://cmake.org/cmake/help/latest/command/find_path.html
  IF (PARPACK_INCLUDE_DIR AND PARPACK_LIB)
    SET(PARPACKLIBS_FOUND TRUE)
    SET(PARPACK_FOUND TRUE)
    SET(PARPACK_LIBRARIES ${PARPACK_LIB})
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
  ENDIF()

ENDIF()


MARK_AS_ADVANCED( ARPACK_LIBRARIES ARPACK_INCLUDE_DIR ARPACK_FOUND ARPACK_LIB
                  ARPACKLIBS_FOUND PARPACK_LIBRARIES PARPACK_INCLUDE_DIR
                  PARPACK_FOUND PARPACK_LIB PARPACKLIBS_FOUND )
