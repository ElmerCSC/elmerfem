
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
  # CMAKE_CURRENT_SOURCE_DIR is this directory in both the in-tree and the
  # standalone build.  The CMAKE_SOURCE_DIR/elmergrid/src spelling that used to
  # be listed here named the same directory when Elmer is the top level project
  # and a nonexistent one otherwise.
  TARGET_INCLUDE_DIRECTORIES(ElmerGrid SYSTEM BEFORE PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/metis-5.1.0/include)
ENDIF()

# config.h is generated into the elmergrid binary directory and reached as
# "../config.h" from the sources here, so this directory is what is needed.
TARGET_INCLUDE_DIRECTORIES(ElmerGrid PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR})

TARGET_COMPILE_DEFINITIONS(ElmerGrid PRIVATE DISABLE_MATC)

TARGET_LINK_LIBRARIES(ElmerGrid PRIVATE metis m)
IF(NOT(WIN32))
  SET_TARGET_PROPERTIES(ElmerGrid PROPERTIES INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING}")
ENDIF()

INSTALL(TARGETS ElmerGrid RUNTIME DESTINATION "bin")

# ---------------------------------------------------------------------------
# ElmerGrid core library, used by the ElmerGrid plugin inside ElmerGUI.
#
# ElmerGUI embeds the ElmerGrid mesh handling routines.  Until now it carried
# its own copies of egconvert.c, egmesh.c, egnative.c and egutils.c under
# ElmerGUI/Application/plugins, renamed to .cpp and refreshed by hand with
# plugins/grid2gui.sh.  Those copies are gone and ElmerGUI links this library
# instead, so there is a single copy of these sources in the tree.
#
# EG_PLUGIN selects the embedded configuration (USE_MATC 0, USE_METIS 0), which
# is what ElmerGUI has always defined for its copies.  It is PUBLIC so that the
# setting reaching this library and the setting reaching ElmerGUI's own sources
# cannot drift apart.
#
# The headers declaring these routines carry extern "C" guards, so C++ callers
# such as ElmerGUI can include them directly.
# ---------------------------------------------------------------------------
SET(elmergrid_core_SRCS
  egconvert.c egmesh.c egnative.c egutils.c
  egconvert.h egmesh.h egnative.h egutils.h
  egdef.h egtypes.h)

ADD_LIBRARY(elmergrid_core STATIC ${elmergrid_core_SRCS})

TARGET_COMPILE_DEFINITIONS(elmergrid_core PUBLIC EG_PLUGIN)

TARGET_INCLUDE_DIRECTORIES(elmergrid_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

SET_TARGET_PROPERTIES(elmergrid_core PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

IF(NOT(WIN32))
  TARGET_LINK_LIBRARIES(elmergrid_core PUBLIC m)
ENDIF()
