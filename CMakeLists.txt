# 
# Experimental CMake configuration script for Elmer
# 
# Authors:
#  Sami Ilvonen, CSC - IT Center for Science, Ltd.
#  Phil Weir, NUMA Engineering Services, Ltd.
#  Mikko Byckling, CSC - IT Center for Science Ltd.
#
# First public beta release 11th of October, 2013
#
# Note that this version is highly experimental and includes
# only a very limited subset of features that are supported
# by the original autotools scripts.
#

PROJECT(Elmer Fortran C CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(APPLE)
  SET(CMAKE_MACOSX_RPATH 1)
ENDIF()

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules")

# Defaults to cache variables
SET(WITH_MPI TRUE CACHE BOOL "Link against real MPI")
SET(WITH_ELMERGUI TRUE CACHE BOOL "Include ElmerGUI")
SET(WITH_ELMERGUITESTER FALSE CACHE BOOL "Include ElmerGUI tester")
SET(ELMER_INSTALL_LIB_DIR "lib/elmersolver" CACHE PATH "Location of elmer shared libraries.")

IF(NOT(CMAKE_BUILD_TYPE MATCHES "Release"))
  ENABLE_TESTING()
ENDIF()

SET(ELMER_FEM_MAJOR_VERSION 7)
SET(ELMER_FEM_MINOR_VERSION 0)
SET(ELMER_FEM_VERSION
  ${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION})

SET(VERSION ${ELMER_FEM_VERSION})
SET(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

IF(WITH_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
ENDIF(WITH_MPI)

MARK_AS_ADVANCED(WITH_TRILINOS)
IF(WITH_TRILINOS)
  FIND_PACKAGE(ML REQUIRED) 
  FIND_PACKAGE(Trilinos REQUIRED)
  FIND_PACKAGE(Belos REQUIRED)
  IF(ML_DIR)
    ADD_DEFINITIONS(-DHAVE_TRILINOS)
  ENDIF()
ENDIF()

# MPI and LAPACK
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

# Mumps
IF(MPI_FOUND)
  FIND_PACKAGE(Mumps)
ENDIF()
IF(MUMPS_FOUND)
  FIND_PACKAGE(SCALAPACK)
  IF(SCALAPACK_FOUND)
     # Add Mumps compilation flags and libraries 
     ADD_DEFINITIONS(-DHAVE_MUMPS)
     INCLUDE_DIRECTORIES(${MUMPS_INCLUDE_PATH})
     SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} ${SCALAPACK_LIBRARIES})
     LINK_DIRECTORIES(${MUMPS_LIBRARIES} ${SCALAPACK_LIBRARIES})
  ELSE()
     MESSAGE(STATUS "SCALAPACK library not found. Mumps not in use.")
  ENDIF()
ENDIF()

# Check if Fortran compiler supports contiguous keyword
INCLUDE(testContiguous)
IF(CMAKE_Fortran_COMPILER_SUPPORTS_CONTIGUOUS)
  ADD_DEFINITIONS(-DCONTIG=,CONTIGUOUS)
ELSE()
  ADD_DEFINITIONS(-DCONTIG=)
ENDIF()

# Check if this is mingw toolchain. 
IF(MINGW)
  ADD_DEFINITIONS(-DMINGW32)
ENDIF()

IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF()

# Check if Fortran compiler supports execute_command_line 
INCLUDE(testExecutecommand)
IF(CMAKE_Fortran_COMPILER_SUPPORTS_EXECUTECOMMANDLINE)
  ADD_DEFINITIONS(-DHAVE_EXECUTECOMMANDLINE)
ENDIF()

ADD_DEFINITIONS(-DUSE_ISO_C_BINDINGS)
SET(FC_STDCALLBULL " ")
SET(ELMER_LINKTYP 1)
#SET(FC_FUNC FC_GLOBAL)
#SET(FC_FUNC_ FC_GLOBAL_)
 
# Add support for arpack 
ADD_DEFINITIONS(-DUSE_ARPACK)

#INCLUDE(FortranCInterface)
#FortranCInterface_HEADER(FCMangle.h MACRO_NAMESPACE "FC_")
#FortranCInterface_VERIFY(CXX)

# Set up rpaths to point ELMER_INSTALL_LIB_DIR
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")   
  SET(ELMERSOLVER_RPATH_STRING "\$ORIGIN/../${ELMER_INSTALL_LIB_DIR}")
ELSEIF(APPLE)                               
  SET(ELMERSOLVER_RPATH_STRING "@loader_path/../${ELMER_INSTALL_LIB_DIR}")
ENDIF()

# Uncomment these to use rpath with linked libraries
# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) 
# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# Remove system directories (from the CMake guide)
# LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
# IF("${isSystemDir}" STREQUAL "-1")
#    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# ENDIF("${isSystemDir}" STREQUAL "-1")
MARK_AS_ADVANCED(USE_EIO)

ADD_SUBDIRECTORY(matc)
IF(WITH_EIO)
  ADD_SUBDIRECTORY(eio)
  SET(HAVE_EIO TRUE)
  MARK_AS_ADVANCED(HAVE_EIO)
ENDIF()
ADD_SUBDIRECTORY(umfpack)
# ADD_SUBDIRECTORY(hutiter)
ADD_SUBDIRECTORY(fhutiter)
ADD_SUBDIRECTORY(meshgen2d)
ADD_SUBDIRECTORY(fem)
ADD_SUBDIRECTORY(mathlibs)
ADD_SUBDIRECTORY(elmergrid)
ADD_SUBDIRECTORY(license_texts)

IF(WITH_ELMERGUI)
  ADD_SUBDIRECTORY(ElmerGUI)
ENDIF()

IF(WITH_ELMERGUITESTER)
  ADD_SUBDIRECTORY(ElmerGUItester)
ENDIF(WITH_ELMERGUITESTER)

#INCLUDE(FeatureSummary)
#FEATURE_SUMMARY(WHAT ENABLED_FEATURES
#  INCLUDE_QUIET_PACKAGES
#  DESCRIPTION "Enabled Features:"
#  VAR enabledFeaturesText)
#MESSAGE(STATUS "${enabledFeaturesText}")

# 


MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  BLAS library:   " ${BLAS_LIBRARIES})
MESSAGE(STATUS "  LAPACK library: " ${LAPACK_LIBRARIES})
MESSAGE(STATUS "------------------------------------------------")
IF(WITH_MPI)
MESSAGE(STATUS "  MPI Fortran:    " ${MPI_Fortran_FOUND})
MESSAGE(STATUS "  MPI Fortran compiler:    " ${MPI_Fortran_COMPILER} 
  ${MPI_Fortran_COMPILER_FLAGS})
MESSAGE(STATUS "  MPI Fortran include dir: " ${MPI_Fortran_INCLUDE_PATH})
MESSAGE(STATUS "  MPI Fortran libraries:   " ${MPI_Fortran_LIBRARIES})
MESSAGE(STATUS "  MPI Fortran link flags:  " ${MPI_Fortran_LINK_FLAGS})
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  MPI C:             " ${MPI_C_FOUND})
MESSAGE(STATUS "  MPI C compiler:    " ${MPI_C_COMPILER} 
  ${MPI_C_COMPILER_FLAGS})
MESSAGE(STATUS "  MPI C include dir: " ${MPI_C_INCLUDE_PATH})
MESSAGE(STATUS "  MPI C libraries:   " ${MPI_C_LIBRARIES})
MESSAGE(STATUS "  MPI C flags:       " ${MPI_C_LINK_FLAGS})
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  Mumps:             " ${MUMPS_FOUND})
MESSAGE(STATUS "  Mumps include:     " ${MUMPS_INCLUDE_PATH})
MESSAGE(STATUS "  Mumps libraries:   " ${MUMPS_LIBRARIES})
MESSAGE(STATUS "------------------------------------------------")
ELSE()
MESSAGE(STATUS "  Fortran compiler:  " ${CMAKE_Fortran_COMPILER})
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  C compiler:        " ${CMAKE_C_COMPILER})
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  C++ compiler:      " ${CMAKE_CXX_COMPILER})
MESSAGE(STATUS "------------------------------------------------")
ENDIF()
IF(CMAKE_BUILD_TYPE STREQUAL "Release")
  MESSAGE(STATUS "  Building release")
ENDIF()
IF(WITH_ELMERGUI)
  MESSAGE(STATUS "  Building ElmerGUI")
ENDIF()
IF(WITH_ELMERGUITESTER)
  MESSAGE(STATUS "  Building ElmerGUI tester")
ENDIF()
MESSAGE(STATUS "------------------------------------------------")

# Packaging
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cpack/ElmerCPack.cmake)
