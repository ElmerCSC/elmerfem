# - cmake script for finding NetCDF libraries

#  ZOLTAN_INCLUDE_DIR  - user modifiable choice of where to zoltan include dir
#  ZOLTAN_LIBRARY    - user modifiable choice of where zoltan library is

# This module returns these variables for the rest of the project to use.
#
#  Zoltan_FOUND              - True if ZOLTAN found
#  ZOLTAN_LIBRARY            - zoltan library is
#  ZOLTAN_INCLUDE_DIR       - zoltan include dir.

IF(ZOLTAN_LIBRARY AND ZOLTAN_INCLUDE_DIR)
  SET(Zoltan_FOUND TRUE)
ELSE()
  SET(Zoltan_FOUND FALSE)
ENDIF()

IF(NOT Zoltan_FOUND)
  MESSAGE(STATUS "Finding Zoltan library")
  # Try to find with CMake config file of upstream Zoltan.
  FIND_PACKAGE(Zoltan CONFIG NAMES Zoltan zoltan)
  IF(Zoltan_FOUND)
    GET_TARGET_PROPERTY(ZOLTAN_INCLUDE_DIR Zoltan::zoltan INTERFACE_INCLUDE_DIRECTORIES)
    GET_TARGET_PROPERTY(ZOLTAN_LIBRARY Zoltan::zoltan IMPORTED_LOCATION_RELEASE)
    # Check if a debug build was used
    IF(NOT ZOLTAN_LIBRARY)
      GET_TARGET_PROPERTY(ZOLTAN_LIBRARY Zoltan::zoltan IMPORTED_LOCATION_DEBUG)
    ENDIF()
    IF(Zoltan_ENABLE_ParMETIS)
      FIND_PACKAGE(ParMetis)
    ENDIF()
  ENDIF()

  IF(NOT ZOLTAN_INCLUDE_DIR AND NOT ZOLTAN_LIBRARY)

    INCLUDE(FindPackageHandleStandardArgs)
    MESSAGE(STATUS "Manual search of Zoltan library")

    SET(ZOLTANINCLUDE
      "${ZOLTANROOT}/include"
      "$ENV{ZOLTANROOT}/include"
      "${ZOLTAN_ROOT}/include"
      "$ENV{ZOLTAN_ROOT}/include"
      "${CMAKE_SOURCE_DIR}/zoltan/include"
      INTERNAL)
    FIND_PATH(ZOLTAN_INCLUDE_DIR zoltan.h HINTS ${ZOLTANINCLUDE})

    SET(ZOLTANLIB
      "${ZOLTANROOT}/lib"
      "$ENV{ZOLTANROOT}/lib"
      "${ZOLTAN_ROOT}/lib"
      "$ENV{ZOLTAN_ROOT}/lib"
      "${CMAKE_SOURCE_DIR}/zoltan/lib"
      INTERNAL)
    FIND_LIBRARY(ZOLTAN_LIBRARY zoltan HINTS ${ZOLTANLIB})

  ENDIF(NOT ZOLTAN_INCLUDE_DIR AND NOT ZOLTAN_LIBRARY)

ENDIF(NOT Zoltan_FOUND)

# This checks could be inadequate because this variables are not empty if nothing found
# Other option is to use the keyword REQUIRED, but this will increase cmake version to 3.18
# https://cmake.org/cmake/help/latest/command/find_library.html
# https://cmake.org/cmake/help/latest/command/find_path.html
IF (ZOLTAN_INCLUDE_DIR AND ZOLTAN_LIBRARY)
  SET(Zoltan_FOUND TRUE)
 #The config script was not used, define the target manually
  IF(NOT TARGET Zoltan::zoltan)
    ADD_LIBRARY(Zoltan::zoltan SHARED IMPORTED)
    SET_TARGET_PROPERTIES(Zoltan::zoltan PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES ${ZOLTAN_INCLUDE_DIR}
                          IMPORTED_LOCATION ${ZOLTAN_LIBRARY} )
  ENDIF()
  IF (NOT Zoltan_FIND_QUIETLY)
    MESSAGE(STATUS "Zoltan library found")
    MESSAGE(STATUS "Zoltan include path: ${ZOLTAN_INCLUDE_DIR}")
    MESSAGE(STATUS "Zoltan library:      ${ZOLTAN_LIBRARY}")
  ENDIF()
ELSE()
  SET(Zoltan_FOUND FALSE)
  IF (Zoltan_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Zoltan library not found.")
  ENDIF()

ENDIF(ZOLTAN_INCLUDE_DIR AND ZOLTAN_LIBRARY)

MARK_AS_ADVANCED( ZOLTANINCLUDE ZOLTANLIB ZOLTAN_INCLUDE_DIR ZOLTAN_LIBRARY )
