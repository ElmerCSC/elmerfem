# CMake script for finding Hypre
# Juhani Kataja, CSC - IT Center for Science Ltd.
# 2014/08
#
# Hint variables:
# * HYPREROOT (env, cmake)
# * HYPRE_ROOT (env, cmake)
# * HYPRE_INCLUDE_DIR (cmake)
# * HYPRE_LIBRARY_DIR (cmake)

cmake_minimum_required(VERSION 2.8)

# If Hypre libraries are already defined, do nothing
IF(Hypre_LIBRARIES)
  IF(Hypre_INCLUDE_DIR)
    SET(Hypre_FOUND TRUE)
    RETURN()
  ENDIF()
ENDIF()

set(Hypre_FOUND FALSE)

find_path(Hypre_INCLUDE_DIR NAMES HYPRE.h
  HINTS
  "${HYPREROOT}/include"
  "$ENV{HYPREROOT}/include"
  "${HYPRE_ROOT}/include"
  "$ENV{HYPRE_ROOT}/include"
  "${HYPRE_INCLUDE_DIR}"
  "${CMAKE_SOURCE_DIR}/hypre/include"
  )

find_library(Hypre_LIBRARY NAMES HYPRE
  HINTS
  "${HYPREROOT}/lib"
  "$ENV{HYPREROOT}/lib"
  "${HYPRE_ROOT}/lib"
  "$ENV{HYPRE_ROOT}/lib"
  "${HYPRE_LIBRARY_DIR}"
  "${CMAKE_SOURCE_DIR}/hypre/lib"
  )

list(APPEND Hypre_LIBRARIES ${Hypre_LIBRARY})

unset(Hypre_LIBRARY CACHE)

foreach(_comp ${Hypre_FIND_COMPONENTS})
  find_library(_Hypre_LIB NAMES HYPRE_${_comp}
    HINTS
    "${HYPREROOT}/lib"
    "$ENV{HYPREROOT}/lib"
  "${HYPRE_LIBRARY_DIR}"
    )
  
  IF(NOT _Hypre_LIB)
    IF(${Hypre_FIND_REQUIRED})
      IF(${Hypre_FIND_REQUIRED_${_comp}})
        MESSAGE(FATAL_ERROR "HYPRE_${_comp}: ${_Hypre_LIB_${_comp}}")
      ENDIF()
    ENDIF()
  ENDIF()

  IF(_Hypre_LIB)
    list(APPEND Hypre_LIBRARIES ${_Hypre_LIB})
  ENDIF()
  unset(_Hypre_LIB CACHE)
endforeach()

IF(Hypre_INCLUDE_DIR AND Hypre_LIBRARIES)
  SET(Hypre_FOUND TRUE)
ENDIF()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Hypre DEFAULT_MSG Hypre_LIBRARIES Hypre_INCLUDE_DIR)

set(Hypre_LIBRARIES ${Hypre_LIBRARIES} CACHE FILEPATH "Hypre Libraries")

MARK_AS_ADVANCED(Hypre_INCLUDE_DIR Hypre_LIBRARIES)
