# CMake script for finding Hypre
# Juhani Kataja, CSC - IT Center for Science Ltd.
# 2014/08

set(Hypre_FOUND FALSE)

find_path(Hypre_INCLUDE_DIR NAMES HYPRE.h
  HINTS
  "${Hypre_DIR}/include"
  "$ENV{HYPREROOT}/include"
  )

find_library(Hypre_LIBRARIES NAMES HYPRE
  HINTS
  "${Hypre_DIR}/lib"
  "$ENV{HYPREROOT}/lib"
  )

IF(${Hypre_LIBRARIES} MATCHES NOTFOUND)
  IF(${Hypre_FIND_REQUIRED})
    MESSAGE(FATAL_ERROR "HYPRE: ${Hypre_LIBRARIES}")
  ELSE()
    IF(NOT Hypre_FIND_QUIETLY)
      MESSAGE(STATUS "Hypre not found.")
    ENDIF()
  ENDIF()
  SET(Hypre_FOUND FALSE)
  RETURN()
ENDIF()

foreach(_comp ${Hypre_FIND_COMPONENTS})
  find_library(_Hypre_LIB NAMES HYPRE_${_comp}
    HINTS
    ${Hypre_DIR})
  
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

IF (Hypre_FOUND)
   IF (NOT Hypre_FIND_QUIETLY)
      MESSAGE(STATUS "A library with Hypre API found.")
      MESSAGE(STATUS "Hypre_INCLUDE_DIR: ${Hypre_INCLUDE_DIR}")
      MESSAGE(STATUS "Hypre_LIBRARIES: ${Hypre_LIBRARIES}")
   ENDIF()
ELSE()
   IF (Hypre_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Hypre libraries not found.")
   ENDIF()
ENDIF()

MARK_AS_ADVANCED(Hypre_INCLUDE_PATH Hypre_LIBRARIES)



