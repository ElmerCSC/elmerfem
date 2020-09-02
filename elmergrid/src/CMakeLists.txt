
IF(CMAKE_CROSSCOMPILING)
  IF(NATIVE_BUILD_FLAGS)
    IF(CMAKE_BUILD_TYPE)
      STRING(TOUPPER _${CMAKE_BUILD_TYPE} _build_type)
      SET(CMAKE_C_FLAGS${_build_type} "")
    ENDIF()
    SET(CMAKE_C_FLAGS ${NATIVE_BUILD_FLAGS})
    MESSAGE(STATUS "Resetting build flags for Elmergrid: ${NATIVE_BUILD_FLAGS}")
  ELSE()
    MESSAGE(WARNING "Crosscompiling ElmerGrid, no native build flags set")
  ENDIF()
ENDIF()

# Don't install metis
SET(METIS_INSTALL FALSE CACHE BOOL "")
ADD_SUBDIRECTORY(metis-5.1.0)

SET(elmergrid_SRCS egdef.h egparallel.h egnative.h egconvert.h
  egexport.h egextra.h egmesh.h 
  egtypes.h egutils.h egparallel.c egnative.c egconvert.c egexport.c
  egextra.c egmesh.c fempre.c 
  egutils.c)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/metis)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/elmergrid/src/metis)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/elmergrid/src)

ADD_DEFINITIONS(-DDISABLE_MATC)

ADD_EXECUTABLE(ElmerGrid ${elmergrid_SRCS})

TARGET_LINK_LIBRARIES(ElmerGrid metis m)
IF(NOT(WIN32))
  SET_TARGET_PROPERTIES(ElmerGrid PROPERTIES INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING}")
ENDIF()

INSTALL(TARGETS ElmerGrid RUNTIME DESTINATION "bin")
