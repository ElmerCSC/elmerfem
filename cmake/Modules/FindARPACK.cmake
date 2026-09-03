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

  # Plain search first; the upstream CMake config package is consulted only if
  # it finds nothing.
  #
  # FIND_PACKAGE(... CONFIG) does not merely look for a config file, it
  # executes it, so a distribution shipping a broken one makes this module fail
  # in a way no test of the result can prevent. Two in circulation today:
  # Debian and Ubuntu install an arpackng-config.cmake whose
  # arpackngTargets.cmake is not packaged, and MSYS2 ships one that declares a
  # parpack target from a separate package that need not be installed. Both
  # abort the configure while reading the file.
  #
  # Searching for the library directly costs nothing when the config package is
  # good, and avoids opening it at all when it is not.
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

  # Nothing found by hand: try the config package after all.
  IF(NOT ARPACK_LIBRARIES OR NOT ARPACK_INCLUDE_DIR)
    # Try to find with CMake config file of upstream arpack.
    FIND_PACKAGE(ARPACK CONFIG NAMES arpack arpackng arpack-ng)
    # IF(TARGET arpack) rather than IF(ARPACK_FOUND): CMake sets ARPACK_FOUND
    # when a config file was located and ran to its end, even if it created no
    # target, and GET_TARGET_PROPERTY on a target that does not exist adds
    # three more errors to whatever the config file already reported. The
    # target is what is actually wanted here.
    IF(TARGET arpack)
      GET_TARGET_PROPERTY(ARPACK_INCLUDE_DIR arpack INTERFACE_INCLUDE_DIRECTORIES)
      GET_TARGET_PROPERTY(ARPACK_LIBRARIES arpack IMPORTED_LOCATION_RELEASE)
      # Check if a debug build type was used
      IF(NOT ARPACK_LIBRARIES)
        GET_TARGET_PROPERTY(ARPACK_LIBRARIES arpack IMPORTED_LOCATION_DEBUG)
      ENDIF()
    ENDIF()
  ENDIF(NOT ARPACK_LIBRARIES OR NOT ARPACK_INCLUDE_DIR)

ENDIF(NOT ARPACK_FOUND)

# This checks could be inadequate because this variables are not empty if nothing found
# Other option is to use the keyword REQUIRED, but this will increase cmake version to 3.18
# https://cmake.org/cmake/help/latest/command/find_library.html
# https://cmake.org/cmake/help/latest/command/find_path.html
IF (ARPACK_INCLUDE_DIR AND ARPACK_LIBRARIES)
  SET(ARPACK_FOUND TRUE)
  #The config script was not used, define the target manually
  #
  # UNKNOWN rather than SHARED. A SHARED imported target needs IMPORTED_IMPLIB
  # on Windows, and only IMPORTED_LOCATION is set here, so anything linking the
  # bare target name got "arpack-NOTFOUND" on the link line while the configure
  # reported ARPACK found. UNKNOWN links IMPORTED_LOCATION directly and is
  # correct whether FIND_LIBRARY returned a .so, a .a or an import library.
  IF(NOT TARGET arpack)
    ADD_LIBRARY(arpack UNKNOWN IMPORTED)
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
