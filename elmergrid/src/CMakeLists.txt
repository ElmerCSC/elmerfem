
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

# ---------------------------------------------------------------------------
# MATC in ElmerGrid
#
# ElmerGrid can evaluate MATC expressions in .grd files: egnative.c calls
# mtc_domath() and mtc_init() under "#if USE_MATC", and egdef.h no longer derives
# USE_MATC from DISABLE_MATC.  Defining DISABLE_MATC unconditionally, as this
# file used to, made that unreachable and left no way to ask for it.
# DISABLE_MATC is no longer used.
#
# The option is OFF by default, so nothing changes unless it is asked for.
# Enabling it needs the matc library, which the Elmer build provides:
# ADD_SUBDIRECTORY(matc) in the top level CMakeLists.txt is unconditional and
# runs before this directory.  The TARGET check keeps a standalone ElmerGrid
# build working, where no such target exists and MATC simply stays off.
#
# No header is needed: egnative.c declares both functions itself.
#
# See ElmerCSC/elmerfem#899 for what enabling this is ultimately for.
# ---------------------------------------------------------------------------
# The flag, WITH_MATC, is input when calling cmake and now controls
# addition of MATC in ElmerGrid, ElmerGUI, and ElmerVTK, in an
# 'all or nothing' fashion.
#
# The flag, USE_MATC_ALL, is set by cmake based on whether MATC
# exists as a cmake target and cmake has been called with WITH_MATC.
#
# The flag, USE_MATC_ALL, is then checked and the appropriate
# target compile definition, USE_MATC, is added for each program.

OPTION(WITH_MATC "Build Elmer with MATC expression support" OFF)
OPTION(USE_MATC_ALL "Compile with MATC in Elmer." OFF)
MARK_AS_ADVANCED(USE_MATC_ALL)

IF(WITH_MATC AND NOT TARGET matc)
  MESSAGE(WARNING
    "  [ElmerGrid] WITH_MATC is ON but no 'matc' target exists; "
    "building without MATC.")
ENDIF()

MESSAGE(STATUS "------------------------------------------------")
IF(WITH_MATC AND TARGET matc)
  MESSAGE(STATUS "  [ElmerGrid] Compiling with MATC")
  SET(USE_MATC_ALL TRUE CACHE BOOL "Enable MATC in Elmer." FORCE)
  TARGET_COMPILE_DEFINITIONS(ElmerGrid PUBLIC USE_MATC=1)
  TARGET_LINK_LIBRARIES(ElmerGrid PRIVATE matc)
ELSE()
  MESSAGE(STATUS "  [ElmerGrid] Compiling without MATC")
  TARGET_COMPILE_DEFINITIONS(ElmerGrid PUBLIC USE_MATC=0)
  SET(USE_MATC_ALL FALSE CACHE BOOL "Disable MATC in Elmer." FORCE)
ENDIF()

TARGET_LINK_LIBRARIES(ElmerGrid PRIVATE metis m)
IF(NOT(WIN32))
  SET_TARGET_PROPERTIES(ElmerGrid PROPERTIES INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING}")
ENDIF()

INSTALL(TARGETS ElmerGrid RUNTIME DESTINATION "bin")

# ---------------------------------------------------------------------------
# ElmerGrid core library, used by the ElmerGrid plugin inside ElmerGUI.
#
# ElmerGUI embeds the ElmerGrid input file and mesh handling routines.  Until now it
# carried its own copies of egconvert.c, egmesh.c, egnative.c and egutils.c under
# ElmerGUI/Application/plugins, renamed to .cpp and refreshed by hand with
# plugins/grid2gui.sh.  Those copies are gone and ElmerGUI links this library
# instead, so there is a single copy of these sources in the tree.
#
# EG_PLUGIN previously selected the embedded configuration (USE_MATC 0, USE_METIS 0),
# which is what ElmerGUI has always defined for its copies.  The flag, EG_PLUGIN, is
# no longer used.  The embedded configuration is now controlled by the flags,
# WITH_MATC and USE_MATC_ALL.
#
# The headers declaring these routines carry extern "C" guards, so C++ callers
# such as ElmerGUI can include them directly.
# ---------------------------------------------------------------------------
SET(elmergrid_core_SRCS
  egconvert.c egmesh.c egnative.c egutils.c
  egconvert.h egmesh.h egnative.h egutils.h
  egdef.h egtypes.h)

ADD_LIBRARY(elmergrid_core STATIC ${elmergrid_core_SRCS})

IF(USE_MATC_ALL)
  MESSAGE(STATUS "  [ElmerGrid_core] Compiling with MATC")
  TARGET_COMPILE_DEFINITIONS(elmergrid_core PUBLIC USE_MATC=1)
ELSE()
  MESSAGE(STATUS "  [ElmerGrid_core] Compiling without MATC")
  TARGET_COMPILE_DEFINITIONS(elmergrid_core PUBLIC USE_MATC=0)
ENDIF()
MESSAGE(STATUS "------------------------------------------------")

TARGET_INCLUDE_DIRECTORIES(elmergrid_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

SET_TARGET_PROPERTIES(elmergrid_core PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

IF(NOT(WIN32))
  TARGET_LINK_LIBRARIES(elmergrid_core PUBLIC m)
ENDIF()
