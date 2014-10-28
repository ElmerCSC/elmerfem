# CMake script for finding ParMetis

# If libraries are already defined, do nothing 
IF(ParMetis_LIBRARIES AND ParMetis_INCLUDE_DIR)
  SET(ParMetis_FOUND TRUE)
  RETURN()
ENDIF()

SET(ParMetis_FOUND FALSE)
MESSAGE(STATUS "Finding ParMetis")

SET(PARMETISINCLUDE 
  "${PARMETISROOT}/include"
  "$ENV{PARMETISROOT}/include"
  "$ENV{PARMETIS_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/parmetis/include"
  INTERNAL)

# Find include
SET(PARMETIS_INCLUDENAME "parmetis.h" INTERNAL)
FIND_PATH(ParMetis_INCLUDE_DIR
  NAMES
  ${PARMETIS_INCLUDENAME} 
  HINTS 
  ${PARMETISINCLUDE}
  )

SET(PARMETISLIB
  "${PARMETISROOT}/lib"
  "$ENV{PARMETISROOT}/lib"
  "$ENV{PARMETIS_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/parmetis/lib"
  INTERNAL)

# Find library
FIND_LIBRARY(ParMetis_LIBRARIES 
  parmetis
  HINTS
  ${PARMETISLIB}
  )

IF (ParMetis_LIBRARIES AND ParMetis_INCLUDE_DIR)
  SET(ParMetis_FOUND TRUE)
ENDIF()

IF (ParMetis_FOUND) 
  IF (NOT ParMetis_FIND_QUIETLY)
    MESSAGE(STATUS "A library with ParMetis API found.")
    MESSAGE(STATUS "ParMetis include dir: ${ParMetis_INCLUDE_DIR}")
    MESSAGE(STATUS "ParMetis libraries: ${ParMetis_LIBRARIES}")
  ENDIF()
ELSE()
  IF (ParMetis_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "ParMetis libraries not found.")
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  PARMETISINCLUDE
  PARMETISLIB
  PARMETIS_INCLUDENAME
  ParMetis_INCLUDE_DIR 
  ParMetis_LIBRARIES 
  )