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
#  YAC_VERSION            - Version of YAC found

# Minimum required YAC version (must have CMake config support)
# Note: Currently a soft requirement with fallback to manual search.
#       Will become a hard requirement (remove fallback) in the future.

SET(YAC_MINIMUM_VERSION "3.12.0")

# # handle the QUIETLY and REQUIRED arguments and set YAC_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

# If YAC is already defined, do nothing
IF(YAC_LIBRARIES AND YAC_INCLUDE_DIR)
   SET(YAC_FOUND TRUE)
   RETURN()
ENDIF()

# First, try to find YAC using its CMake config files (for modern CMake-based installations)
SET(YAC_CONFIG_PATHS
  "${YACROOT}"
  "$ENV{YACROOT}"
  "${YAC_ROOT}"
  "$ENV{YAC_ROOT}"
)

# Provide glue code for YAC's config mode to find dependencies
# YAC's CMake config expects specific variable names

# NetCDF mapping
IF(NETCDF_LIBRARY AND NOT NetCDF_C_LIBRARIES)
  SET(NetCDF_C_LIBRARIES "${NETCDF_LIBRARY}")
ENDIF()
IF(NETCDFF_LIBRARY AND NOT NetCDF_Fortran_LIBRARIES)
  SET(NetCDF_Fortran_LIBRARIES "${NETCDFF_LIBRARY}")
ENDIF()
IF(NETCDF_INCLUDE_DIR AND NOT NetCDF_INCLUDE_DIR)
  SET(NetCDF_INCLUDE_DIR "${NETCDF_INCLUDE_DIR}")
ENDIF()

# YAXT mapping
IF(YAXT_LIBRARY AND NOT YAXT::YAXT_C)
  # Create imported target for YAXT if it doesn't exist
  IF(NOT TARGET YAXT::YAXT_C)
    ADD_LIBRARY(YAXT::YAXT_C UNKNOWN IMPORTED)
    SET_TARGET_PROPERTIES(YAXT::YAXT_C PROPERTIES
      IMPORTED_LOCATION "${YAXT_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${YAXT_INCLUDE_DIR}")
  ENDIF()
ENDIF()

# libfyaml (YAML) mapping
IF(YAML_LIBRARY AND NOT libfyaml_LIBRARIES)
  SET(libfyaml_LIBRARIES "${YAML_LIBRARY}")
ENDIF()
IF(YAML_INCLUDE_DIR AND NOT libfyaml_INCLUDE_DIRS)
  SET(libfyaml_INCLUDE_DIRS "${YAML_INCLUDE_DIR}")
ENDIF()
IF(YAML_FOUND AND NOT libfyaml_FOUND)
  SET(libfyaml_FOUND TRUE)
ENDIF()

FIND_PACKAGE(YAC CONFIG QUIET HINTS ${YAC_CONFIG_PATHS})

IF(YAC_FOUND)
  # Check YAC version if available
  IF(YAC_VERSION)
    IF(YAC_VERSION VERSION_LESS YAC_MINIMUM_VERSION)
      message(WARNING "Found YAC version ${YAC_VERSION}, but minimum required is ${YAC_MINIMUM_VERSION}")
      SET(YAC_FOUND FALSE)
      # RETURN()  # TODO turn this on when YAC 3.15 with CMake has been released.
    ENDIF()
  ENDIF()

  # YAC was found via config mode - collect all YAC component targets
  SET(YAC_COMPONENT_TARGETS "")

  # YAC exports multiple component targets in the yac:: namespace
  FOREACH(component yac_core yac_mci yac_utils yac_mtime)
    IF(TARGET yac::${component})
      LIST(APPEND YAC_COMPONENT_TARGETS yac::${component})
    ENDIF()
  ENDFOREACH()

  IF(YAC_COMPONENT_TARGETS)
    # For modern CMake, set YAC_LIBRARIES to all component targets
    # This allows: target_link_libraries(... ${YAC_LIBRARIES}) to work correctly
    SET(YAC_LIBRARIES ${YAC_COMPONENT_TARGETS})

    # Get include directories from the first component target
    LIST(GET YAC_COMPONENT_TARGETS 0 FIRST_TARGET)
    GET_TARGET_PROPERTY(YAC_INCLUDE_DIR ${FIRST_TARGET} INTERFACE_INCLUDE_DIRECTORIES)

    IF(YAC_VERSION)
      message(STATUS "Found YAC ${YAC_VERSION} via config mode: targets ${YAC_LIBRARIES}")
    ELSE()
      message(STATUS "Found YAC via config mode: targets ${YAC_LIBRARIES}")
    ENDIF()
    RETURN()
  ELSE()
    message(WARNING "FindYAC: YAC config found but no component targets available")
  ENDIF()
ENDIF()

message(WARNING "Did not find YAC via config mode. "
  "YAC >= ${YAC_MINIMUM_VERSION} with CMake config support is required."
  "Falling back to manual search for older YAC installations without CMake.")

message(WARNING "Fallback to YAC without CMake config support will be removed "
  "in the future, so please upgrade YAC to a modern version with CMake support.")

# =============================================================================
# FALLBACK: Manual search for older YAC installations without CMake support
# =============================================================================
# TODO: Remove this fallback section when making YAC >= 3.15 a hard requirement.
#       YAC 3.15+ provides CMake config files and should be found via CONFIG mode above.
#       The manual search below is only for backwards compatibility with older installations.

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
