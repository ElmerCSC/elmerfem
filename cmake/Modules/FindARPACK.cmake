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

  # The config package is preferred and the plain search below is only the
  # fallback, because an imported target carries more than a file name: usage
  # requirements, interface definitions, and the transitive link dependencies
  # a static build needs. None of that can be reconstructed from a path
  # returned by FIND_LIBRARY.
  #
  # The difficulty is that FIND_PACKAGE(... CONFIG) does not look for a config
  # file, it executes one, and two distributions ship a config that raises
  # FATAL_ERROR when executed. Neither can be caught in this process, and they
  # fail for different reasons, which is why checking that files exist is not
  # enough to tell a broken installation from a sound one:
  #
  #   Debian, Ubuntu -- arpackng-config.cmake include()s an
  #     arpackngTargets.cmake that is not packaged, so the include() fails.
  #
  #   MSYS2 -- every config file is present, but arpackngTargets.cmake declares
  #     an imported parpack target whose library ships in a separate package.
  #     Its own _cmake_import_check_files_for_parpack loop then raises
  #     FATAL_ERROR when only arpack is installed.
  #
  # So the config is exercised in a throwaway CMake process first, and loaded
  # here only if that process survives it. Where the installation is sound,
  # which is the common case, the imported target is used exactly as before.
  #
  # Diagnosed by Juha Ruokolainen on Ubuntu in the discussion on PR #844, and
  # reproduced on MSYS2. The probe runs once and its verdict is cached.
  IF(NOT DEFINED ARPACK_CONFIG_USABLE)
    SET(_arpack_probe "${CMAKE_BINARY_DIR}/CMakeFiles/arpack_config_probe")
    FILE(WRITE "${_arpack_probe}/CMakeLists.txt"
      "cmake_minimum_required(VERSION 3.13)\n"
      "project(arpack_config_probe C)\n"
      "find_package(ARPACK CONFIG NAMES arpack arpackng arpack-ng REQUIRED)\n"
      "if(NOT TARGET arpack)\n"
      "  message(FATAL_ERROR \"config package produced no arpack target\")\n"
      "endif()\n")
    # Same generator, compiler and search paths, or the probe answers a
    # question about a different configuration than the one being configured.
    EXECUTE_PROCESS(
      COMMAND "${CMAKE_COMMAND}"
              -S "${_arpack_probe}"
              -B "${_arpack_probe}/build"
              -G "${CMAKE_GENERATOR}"
              "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
              "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
              "-DARPACK_DIR=${ARPACK_DIR}"
      RESULT_VARIABLE _arpack_probe_rc
      OUTPUT_VARIABLE _arpack_probe_out
      ERROR_VARIABLE  _arpack_probe_out)
    IF(_arpack_probe_rc EQUAL 0)
      SET(ARPACK_CONFIG_USABLE TRUE CACHE INTERNAL "arpack config package loads without error")
    ELSE()
      SET(ARPACK_CONFIG_USABLE FALSE CACHE INTERNAL "arpack config package loads without error")
      MESSAGE(STATUS
        "The installed arpack CMake config package could not be loaded; falling "
        "back to a plain library search. Anything the config would have carried "
        "beyond the library path is lost.")
    ENDIF()
  ENDIF()

  IF(ARPACK_CONFIG_USABLE)
    # Try to find with CMake config file of upstream arpack.
    FIND_PACKAGE(ARPACK CONFIG NAMES arpack arpackng arpack-ng)
  ENDIF()

  # IF(TARGET arpack) rather than IF(ARPACK_FOUND): CMake sets ARPACK_FOUND
  # when a config file was located and ran to its end, even if it created no
  # target, and GET_TARGET_PROPERTY on a target that does not exist adds three
  # more errors to whatever the config file already reported. When the target
  # does exist it is what gets linked; these variables are for reporting and
  # for CMAKE_REQUIRED_LIBRARIES.
  IF(TARGET arpack)
    GET_TARGET_PROPERTY(ARPACK_INCLUDE_DIR arpack INTERFACE_INCLUDE_DIRECTORIES)
    GET_TARGET_PROPERTY(ARPACK_LIBRARIES arpack IMPORTED_LOCATION_RELEASE)
    # Check if a debug build type was used
    IF(NOT ARPACK_LIBRARIES)
      GET_TARGET_PROPERTY(ARPACK_LIBRARIES arpack IMPORTED_LOCATION_DEBUG)
    ENDIF()
    IF(NOT ARPACK_LIBRARIES)
      GET_TARGET_PROPERTY(ARPACK_LIBRARIES arpack IMPORTED_LOCATION)
    ENDIF()
  ENDIF()

  # No usable config package: search for the library by hand.
  IF(NOT ARPACK_LIBRARIES OR NOT ARPACK_INCLUDE_DIR)
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
