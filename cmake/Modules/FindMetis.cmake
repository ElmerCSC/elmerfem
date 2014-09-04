# CMake script for finding Metis

# If libraries are already defined, do nothing 
IF(Metis_LIBRARIES AND Metis_INCLUDE_DIR)
  SET(Metis_FOUND TRUE)
  RETURN()
ENDIF()

SET(Metis_FOUND FALSE)
MESSAGE(STATUS "Finding Metis")

# Try to find Metis
SET(METISINCLUDE 
  "${METISROOT}/include"
  "$ENV{METISROOT}/include"
  "$ENV{METIS_ROOT}/include"
  "$ENV{PARMETISROOT}/include"
  "$ENV{PARMETIS_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/metis/include"
  INTERNAL)

SET(METIS_INCLUDENAME "metis.h" "parmetis.h" INTERNAL)
FIND_PATH(Metis_INCLUDE_DIR
  NAMES
  ${METIS_INCLUDENAME} 
  HINTS 
  ${METISINCLUDE}
  )

SET(METISLIB
  "${METISROOT}/lib"
  "$ENV{METISROOT}/lib"
  "$ENV{METIS_ROOT}/lib"
  "$ENV{PARMETISROOT}/lib"
  "$ENV{PARMETIS_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/metis/lib"
  INTERNAL)

FIND_LIBRARY(Metis_LIBRARIES 
  metis
  HINTS
  ${METISLIB})

IF (Metis_LIBRARIES AND Metis_INCLUDE_DIR)
  SET(Metis_FOUND TRUE)
ENDIF()

IF (Metis_FOUND) 
  IF (NOT Metis_FIND_QUIETLY)
    MESSAGE(STATUS "A library with Metis API found.")
    MESSAGE(STATUS "Metis include dir: ${Metis_INCLUDE_DIR}")
    MESSAGE(STATUS "Metis libraries: ${Metis_LIBRARIES}")
  ENDIF()
ELSE()
  IF (Metis_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Metis libraries not found.")
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  METISINCLUDE
  METISLIB
  METIS_INCLUDENAME
  Metis_INCLUDE_DIR 
  Metis_LIBRARIES 
  )