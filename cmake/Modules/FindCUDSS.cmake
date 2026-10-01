# CMake script for finding NVIDIA's cuDSS GPU sparse direct solver library.
#
# Prefers cuDSS's own cudss-config.cmake, because the CUDA network repo's
# RPM/DEB packages namespace the install path by CUDA major version so
# cuda-12 and cuda-13 can coexist (/usr/lib64/libcudss/12/libcudss.so), not
# the flat <prefix>/include, <prefix>/lib layout FIND_PATH / FIND_LIBRARY
# expects. Falls back to a flat search, for NVIDIA's plain tar.xz archives.
#
# This module returns:
#
#  CUDSS_FOUND          - True if cuDSS was found
#  CUDSS_INCLUDE_DIR    - directory containing cudss.h
#  CUDSS_LIBRARIES      - pass to TARGET_LINK_LIBRARIES: the "cudss" imported
#                          target (config path) or a library path (fallback)
#
# Search hint:
#   - RPM/DEB: the glob below finds /usr/lib64/libcudss/*/cmake/cudss, taking
#     the first match if several CUDA versions are installed; pass
#     -Dcudss_DIR=<...>/cmake/cudss to choose one.

# If cuDSS libraries are already defined, do nothing
IF(CUDSS_LIBRARIES AND CUDSS_INCLUDE_DIR)
   SET(CUDSS_FOUND TRUE)
   RETURN()
ENDIF()

INCLUDE(FindPackageHandleStandardArgs)

# ── Preferred: NVIDIA's own package config ──────────────────────────────────
IF(NOT cudss_DIR)
  FILE(GLOB _cudss_config_hints
    "${CUDSS_ROOT}/lib64/libcudss/*/cmake/cudss"
    "${CUDSS_ROOT}/lib/libcudss/*/cmake/cudss"
    "$ENV{CUDSS_ROOT}/lib64/libcudss/*/cmake/cudss"
    "$ENV{CUDSS_ROOT}/lib/libcudss/*/cmake/cudss"
    "/usr/lib64/libcudss/*/cmake/cudss"
    "/usr/lib/libcudss/*/cmake/cudss"
    "/usr/lib/*/libcudss/*/cmake/cudss")
  IF(_cudss_config_hints)
    LIST(GET _cudss_config_hints 0 cudss_DIR)
  ENDIF()
  UNSET(_cudss_config_hints)
ENDIF()

FIND_PACKAGE(cudss CONFIG QUIET)

IF(cudss_FOUND)
  SET(CUDSS_FOUND TRUE)
  SET(CUDSS_INCLUDE_DIR ${cudss_INCLUDE_DIR})
  SET(CUDSS_LIBRARIES cudss) # the imported target; TARGET_LINK_LIBRARIES takes it directly
  MARK_AS_ADVANCED(cudss_DIR)
  IF(NOT CUDSS_FIND_QUIETLY)
    MESSAGE(STATUS "Found cuDSS via upstream cudss-config.cmake:")
    MESSAGE(STATUS "  Version:     ${cudss_VERSION}")
    MESSAGE(STATUS "  Include dir: ${cudss_INCLUDE_DIR}")
    MESSAGE(STATUS "  Library dir: ${cudss_LIBRARY_DIR}")
  ENDIF()
  RETURN()
ENDIF()

# ── Fallback: manual search (flat tarball install) ──────────────────────────
# Not tested, the Elmer GPU containers use the RPM/DEB packages.
SET(CUDSS_FOUND FALSE)
MESSAGE(STATUS "cuDSS package config not found; falling back to a flat-layout search")

SET(CUDSSINCLUDE
  "${CUDSSROOT}/include"
  "$ENV{CUDSSROOT}/include"
  "${CUDSS_ROOT}/include"
  "$ENV{CUDSS_ROOT}/include"
  "${CUDAToolkit_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/cudss/include"
  INTERNAL)

FIND_PATH(CUDSS_INCLUDE_DIR
  NAMES cudss.h
  HINTS ${CUDSSINCLUDE})

SET(CUDSSLIB
  "${CUDSSROOT}/lib"
  "$ENV{CUDSSROOT}/lib"
  "${CUDSSROOT}/lib64"
  "$ENV{CUDSSROOT}/lib64"
  "${CUDSS_ROOT}/lib"
  "$ENV{CUDSS_ROOT}/lib"
  "${CUDSS_ROOT}/lib64"
  "$ENV{CUDSS_ROOT}/lib64"
  "${CUDAToolkit_ROOT}/lib64"
  "${CMAKE_SOURCE_DIR}/cudss/lib"
  INTERNAL)

FIND_LIBRARY(CUDSS_LIBRARY NAMES cudss HINTS ${CUDSSLIB})

IF(CUDSS_INCLUDE_DIR AND CUDSS_LIBRARY)
  UNSET(CUDSS_FAILMSG)
  SET(CUDSS_LIBRARIES ${CUDSS_LIBRARY})
ELSE()
  SET(CUDSS_FAILMSG
    "cuDSS not found. For an RPM/DEB (CUDA network repo) install, pass "
    "-Dcudss_DIR=<...>/cmake/cudss. For a tarball install, set CUDSS_ROOT "
    "to the extracted archive's top-level directory "
    "(containing include/cudss.h and lib/libcudss.so).")
ENDIF()

IF(NOT CUDSS_FAILMSG)
  SET(CUDSS_FOUND TRUE)
ENDIF()

IF(CUDSS_FOUND)
  IF(NOT CUDSS_FIND_QUIETLY)
    MESSAGE(STATUS "Found cuDSS:")
    MESSAGE(STATUS "  Include dir: ${CUDSS_INCLUDE_DIR}")
    MESSAGE(STATUS "  Libraries:   ${CUDSS_LIBRARIES}")
  ENDIF()
ELSE()
  IF(CUDSS_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${CUDSS_FAILMSG})
  ELSE()
    MESSAGE(STATUS "cuDSS not found: ${CUDSS_FAILMSG}")
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  CUDSSINCLUDE CUDSSLIB
  CUDSS_FAILMSG
  CUDSS_INCLUDE_DIR CUDSS_LIBRARY CUDSS_LIBRARIES)
