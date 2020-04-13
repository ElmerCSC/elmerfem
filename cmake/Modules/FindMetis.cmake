# CMake script for finding Metis

# If libraries are already defined, do nothing
if(Metis_LIBRARIES AND Metis_INCLUDE_DIR)
  set(Metis_FOUND TRUE)
  return()
endif()

set(Metis_FOUND FALSE)
message(STATUS "Finding Metis")

# Try to find Metis
set(METISINCLUDE
    "${METISROOT}/include"
    "$ENV{METISROOT}/include"
    "$ENV{METIS_ROOT}/include"
    "$ENV{PARMETISROOT}/include"
    "$ENV{PARMETIS_ROOT}/include"
    "${CMAKE_SOURCE_DIR}/metis/include"
    INTERNAL)

set(METIS_INCLUDENAME "metis.h" "parmetis.h" INTERNAL)
find_path(
  Metis_INCLUDE_DIR metis.h
  PATHS /usr/include
  NAMES ${METIS_INCLUDENAME}
  HINTS ${METISINCLUDE})

set(METISLIB
    "${METISROOT}/lib"
    "$ENV{METISROOT}/lib"
    "$ENV{METIS_ROOT}/lib"
    "$ENV{PARMETISROOT}/lib"
    "$ENV{PARMETIS_ROOT}/lib"
    "${CMAKE_SOURCE_DIR}/metis/lib"
    INTERNAL)

find_library(Metis_LIBRARIES metis HINTS ${METISLIB})

if(Metis_LIBRARIES AND Metis_INCLUDE_DIR)
  set(Metis_FOUND TRUE)
endif()

if(Metis_FOUND)
  if(NOT Metis_FIND_QUIETLY)
    message(STATUS "A library with Metis API found.")
    message(STATUS "Metis include dir: ${Metis_INCLUDE_DIR}")
    message(STATUS "Metis libraries: ${Metis_LIBRARIES}")
  endif()
else()
  if(Metis_FIND_REQUIRED)
    message(FATAL_ERROR "Metis libraries not found.")
  endif()
endif()

mark_as_advanced(METISINCLUDE METISLIB METIS_INCLUDENAME Metis_INCLUDE_DIR
                 Metis_LIBRARIES)
