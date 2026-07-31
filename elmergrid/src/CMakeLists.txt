
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

SET(elmergrid_SRCS egdef.h egparallel.h egnative.h egconvert.h
  egexport.h egextra.h egmesh.h 
  egtypes.h egutils.h egparallel.c egnative.c egconvert.c egexport.c
  egextra.c egmesh.c fempre.c 
  egutils.c)

ADD_EXECUTABLE(ElmerGrid ${elmergrid_SRCS})

IF(EXTERNAL_METIS)
  # use Metis that isinstalled on the system
  FIND_PACKAGE(Metis)
  IF(NOT Metis_FOUND)
    MESSAGE(FATAL_ERROR "  Metis library not found. ")
  ENDIF()
ELSE()
  # do not install bundled metis
  SET(METIS_INSTALL FALSE CACHE BOOL "")
  ADD_SUBDIRECTORY(metis-5.1.0)

  # add metis as system include directory
  TARGET_INCLUDE_DIRECTORIES(ElmerGrid SYSTEM BEFORE PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/metis-5.1.0/include
    ${CMAKE_SOURCE_DIR}/elmergrid/src/metis-5.1.0/include)
ENDIF()

TARGET_INCLUDE_DIRECTORIES(ElmerGrid PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/elmergrid/src)

TARGET_COMPILE_DEFINITIONS(ElmerGrid PRIVATE DISABLE_MATC)

TARGET_LINK_LIBRARIES(ElmerGrid PRIVATE metis m)
IF(NOT(WIN32))
  SET_TARGET_PROPERTIES(ElmerGrid PROPERTIES INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING}")
ENDIF()

INSTALL(TARGETS ElmerGrid RUNTIME DESTINATION "bin")
