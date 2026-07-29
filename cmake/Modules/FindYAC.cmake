# - cmake script for finding YAC libraries
#
# YAC can be installed either via autotools (classic layout: static
# libraries libyac.a, libyac_core.a, libyac_mci.a, libyac_mtime.a next to
# each other) or via CMake (which additionally installs a YAC package
# config, e.g. <prefix>/lib/cmake/yac/yac-config.cmake, exporting imported
# targets such as YAC::yac_mci).
#
# This module first tries YAC's own CMake package config in CONFIG mode
# (this cannot recurse back into this module, since MODULE mode is
# excluded). That correctly captures YAC's internal dependency graph
# (yac_core/yac_pak/yac_utils/yac_mtime/libfyaml/MPI/...) via imported
# targets instead of guessing static library names. If no such config is
# found (classic autotools install), it falls back to the manual search
# below, unchanged.

#  YAC_INCLUDE_DIR  - user modifiable choice of where yac headers are
#  YAC_LIBRARY      - user modifiable choice of where yac library if
#  or YAC_LIBRARIES - user modifiable choice of where yac libraries is

# his module returns these variables for the rest of the project to use.
#
#  YAC_FOUND              - True if YAC found including required interfaces (see below)
#  YAC_LIBRARIES          - All yac related libraries.
#  YAC_INCLUDE_DIR        - All directories to include.

# # handle the QUIETLY and REQUIRED arguments and set YAC_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If YAC is already defined, do nothing
IF(YAC_LIBRARIES AND YAC_INCLUDE_DIR)
   SET(YAC_FOUND TRUE)
   RETURN()
ENDIF()

SET(YAC_FOUND FALSE)

# ------------------------------------------------------------------------
# 1) YAC built with CMake: use its own package config.
# ------------------------------------------------------------------------
FIND_PACKAGE(YAC CONFIG QUIET)

IF (YAC_FOUND AND TARGET YAC::yac_mci)
  GET_TARGET_PROPERTY(YAC_INCLUDE_DIR YAC::yac_mci INTERFACE_INCLUDE_DIRECTORIES)
  SET(YAC_LIBRARIES YAC::yac_mci)
  MARK_AS_ADVANCED(YAC_INCLUDE_DIR YAC_LIBRARIES)
  RETURN()
ENDIF()

# Config mode either found no YAC package config at all (classic autotools
# install, nothing shipped under lib/cmake/yac) or found one without the
# Fortran interface target; fall back to the manual search either way.
SET(YAC_FOUND FALSE)

# ------------------------------------------------------------------------
# 2) YAC built with autotools: manual search for the classic static
#    libraries.
# ------------------------------------------------------------------------
SET(YACINCLUDE
  "${YACROOT}/include"
  "$ENV{YACROOT}/include"
  "${YAC_ROOT}/include"
  "$ENV{YAC_ROOT}/include"
  INTERNAL
  )

FIND_PATH(YAC_INCLUDE_DIR
  yac.h
  HINTS
  ${YACINCLUDE}
  )

SET(YACLIB
  "${YACROOT}/lib"
  "$ENV{YACROOT}/lib"
  "${YAC_ROOT}/lib"
  "$ENV{YAC_ROOT}/lib"
  INTERNAL)

FIND_LIBRARY(YAC_LIBRARY yac HINTS ${YACLIB})

IF (YAC_INCLUDE_DIR AND YAC_LIBRARY)
  UNSET(YAC_FAILMSG)
  SET(YAC_INCLUDE_DIR ${YAC_INCLUDE_DIR})

  # List of YAC static libraries
  SET(YACLIB_NAMES
    libyac.a
    libyac_core.a
    libyac_mci.a
    libyac_mtime.a
  )

  # Determine directory of the found YAC library
  get_filename_component(YACLIB_DIR "${YAC_LIBRARY}" DIRECTORY)

  SET(YAC_LIBRARIES "")
  FOREACH(LIB ${YACLIB_NAMES})
    SET(YAC_LIB_PATH "${YACLIB_DIR}/${LIB}")
    IF (EXISTS "${YAC_LIB_PATH}")
      LIST(APPEND YAC_LIBRARIES "${YAC_LIB_PATH}")
    ELSE()
      # if any of the expected libraries is not found return with error
      SET(YAC_FAILMSG "YAC:           Expected library not found: ${YAC_LIB_PATH}.")
      SET(YACLIB_FOUND FALSE)
      RETURN()
    ENDIF()
  ENDFOREACH()
  SET(YACLIB_FOUND TRUE)
ELSE()
  SET(YAC_FAILMSG "YAC libraries not found.")
ENDIF()

IF (NOT YAC_FAILMSG)
  SET(YAC_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  YACINCLUDE
  YACLIB
  YAC_FAILMSG
  YAC_INCLUDE_DIR
  YAC_LIBRARIES
  YAC_LIBRARY
)
