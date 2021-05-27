# cmake script for finding MUMPS sparse direct solver
INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# If Mumps libraries are already defined, do nothing
IF(Mumps_LIBRARIES AND Mumps_INCLUDE_DIR)
   SET(Mumps_FOUND TRUE)
   RETURN()
ENDIF()

SET(Mumps_FOUND FALSE)
MESSAGE(STATUS "Finding Mumps")

SET(MUMPSINCLUDE
  "${MUMPSROOT}/include"
  "$ENV{MUMPSROOT}/include"
  "${MUMPS_ROOT}/include"
  "$ENV{MUMPS_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/mumps/include"
  INTERNAL
  )
# Try to find Mumps
FIND_PATH(Mumps_INCLUDE_DIR 
  dmumps_struc.h 
  HINTS 
  ${MUMPSINCLUDE}
  )

SET(MUMPSLIB 
  "${MUMPSROOT}/lib"
  "$ENV{MUMPSROOT}/lib"
  "${MUMPS_ROOT}/lib"
  "$ENV{MUMPS_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/mumps/lib"
  INTERNAL)

FIND_LIBRARY(MUMPS_D_LIB dmumps HINTS ${MUMPSLIB})
FIND_LIBRARY(MUMPS_COMMON_LIB mumps_common HINTS ${MUMPSLIB})
FIND_LIBRARY(MUMPS_PORD_LIB pord HINTS ${MUMPSLIB})

IF (Mumps_INCLUDE_DIR AND MUMPS_D_LIB AND MUMPS_COMMON_LIB AND MUMPS_PORD_LIB)
  UNSET(MUMPS_FAILMSG)
  SET(MUMPSLIBS_FOUND TRUE)
  SET(Mumps_LIBRARIES ${MUMPS_D_LIB} ${MUMPS_COMMON_LIB} ${MUMPS_PORD_LIB})
ELSE()
  SET(MUMPS_FAILMSG "Mumps library not found.")
ENDIF()
   
IF (MUMPSLIBS_FOUND)
  # Parallel Mumps always needs Scalapack
  FIND_PACKAGE(SCALAPACK QUIET)
  
  IF(SCALAPACK_FOUND)    
    # Add Mumps compilation flags and libraries 
    SET(Mumps_LIBRARIES ${Mumps_LIBRARIES} ${SCALAPACK_LIBRARIES})

    MESSAGE(STATUS "Checking if Metis library is needed by Mumps")
    UNSET(METIS_NODEND_OUTPUT)
    UNSET(METIS_NODEND_ERROR)
    # Check for Metis
    EXECUTE_PROCESS(COMMAND ${CMAKE_NM} ${MUMPS_D_LIB}
      OUTPUT_VARIABLE METIS_NODEND_OUTPUT
      ERROR_VARIABLE METIS_NODEND_ERROR)
    STRING(FIND "${METIS_NODEND_OUTPUT}" "metis_nodend" METISREF_FOUND)
    
    IF("${METISREF_FOUND}" STREQUAL "-1" AND
	"${METIS_NODEND_ERROR}" STREQUAL "")
      SET(METIS_NEEDED FALSE)
    ELSE()
      SET(METIS_NEEDED TRUE)
    ENDIF()
    
    IF(METIS_NEEDED)
      MESSAGE(STATUS "Checking if Metis library is needed by Mumps -- yes")
      FIND_PACKAGE(Metis QUIET)
    ELSE()
      MESSAGE(STATUS "Checking if Metis library is needed by Mumps -- no")
    ENDIF()
    
    IF(METIS_NEEDED)
      IF(Metis_FOUND)
	SET(Mumps_LIBRARIES ${Mumps_LIBRARIES} ${Metis_LIBRARIES})
	SET(Mumps_INCLUDE_DIR ${Mumps_INCLUDE_DIR} ${Metis_INCLUDE_DIR})
      
	LIST(REMOVE_DUPLICATES Mumps_LIBRARIES)
	LIST(REMOVE_DUPLICATES Mumps_INCLUDE_DIR)
      ELSE()
	SET(MUMPS_FAILMSG 
	  "Metis library not found, needed by found Mumps library.")
      ENDIF()
    ENDIF()

    # Check for ParMetis 
    MESSAGE(STATUS "Checking if ParMetis library is needed by Mumps")
    UNSET(PARMETIS_NODEND_OUTPUT)
    UNSET(PARMETIS_NODEND_ERROR)
    EXECUTE_PROCESS(COMMAND ${CMAKE_NM} ${MUMPS_COMMON_LIB}
      OUTPUT_VARIABLE PARMETIS_NODEND_OUTPUT
      ERROR_VARIABLE PARMETIS_NODEND_ERROR)
    STRING(FIND "${PARMETIS_NODEND_OUTPUT}" "ParMETIS_V3_NodeND" PARMETISREF_FOUND)
    
    IF("${PARMETISREF_FOUND}" STREQUAL "-1" AND 
	"${PARMETIS_NODEND_ERROR}" STREQUAL "")
      SET(PARMETIS_NEEDED FALSE)
    ELSE()
      SET(PARMETIS_NEEDED TRUE)
    ENDIF()
    
    IF(PARMETIS_NEEDED)
      MESSAGE(STATUS "Checking if ParMetis library is needed by Mumps -- yes")
      FIND_PACKAGE(ParMetis QUIET)
    ELSE()
      MESSAGE(STATUS "Checking if ParMetis library is needed by Mumps -- no")
    ENDIF()
    # TODO: Check for PT-Scotch

    IF(PARMETIS_NEEDED)
      IF(ParMetis_FOUND)
	SET(Mumps_LIBRARIES ${Mumps_LIBRARIES} ${ParMetis_LIBRARIES})
	SET(Mumps_INCLUDE_DIR ${Mumps_INCLUDE_DIR} ${ParMetis_INCLUDE_DIR})
	
	LIST(REMOVE_DUPLICATES Mumps_LIBRARIES)
	LIST(REMOVE_DUPLICATES Mumps_INCLUDE_DIR)
      ELSE()
	SET(MUMPS_FAILMSG 
	  "ParMetis library not found, needed by found Mumps library.")
      ENDIF()
    ENDIF()
  
  ELSE()
    SET(MUMPS_FAILMSG 
      "SCALAPACK library not found, required by Mumps.")
  ENDIF()
ENDIF()
 
IF (NOT MUMPS_FAILMSG)
  SET(Mumps_FOUND TRUE)
ENDIF()

IF (Mumps_FOUND)
  IF (NOT Mumps_FIND_QUIETLY)
    MESSAGE(STATUS "A library with Mumps API found.")
    MESSAGE(STATUS "Mumps include dir: ${Mumps_INCLUDE_DIR}")
    MESSAGE(STATUS "Mumps libraries: ${Mumps_LIBRARIES}")
  ENDIF()
ELSE()
  IF (Mumps_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${MUMPS_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  MUMPSINCLUDE
  MUMPSLIB
  MUMPS_FAILMSG
  Mumps_INCLUDE_DIR 
  Mumps_LIBRARIES 
  MUMPS_COMMON_LIB
  MUMPS_D_LIB 
  MUMPS_PORD_LIB 
  SCALAPACK_LIBRARIES)