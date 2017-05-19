# 
# CMake configuration script for Elmer
# 
# Authors:
#  Sami Ilvonen, CSC - IT Center for Science, Ltd.
#  Phil Weir, NUMA Engineering Services, Ltd.
#  Mikko Byckling, CSC - IT Center for Science Ltd.
#  Juhani Kataja, CSC - IT Center for Science Ltd.
#
# First public beta release 11th of October, 2013
#

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type (Release, Debug, RelWithDebugInfo, MinSizeRel)")
ENDIF()

PROJECT(Elmer Fortran C CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)

IF(APPLE)
  SET(CMAKE_MACOSX_RPATH 1) 
ENDIF()

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" ${CMAKE_ROOT}/Modules)

# Defaults to cache variables
SET(WITH_MPI TRUE CACHE BOOL "Use MPI parallelization")
SET(MPI_TEST_MAXPROC 8 CACHE STRING "Maximum number of tasks used in parallel tests")
SET(MPI_TEST_MINPROC 1 CACHE STRING "Minimum number of tasks used in parallel tests")
SET(WITH_OpenMP FALSE CACHE BOOL "Use OpenMP thread parallelization")
SET(WITH_MKL FALSE CACHE BOOL "Use Intel Math Kernel library")
SET(WITH_Mumps FALSE CACHE BOOL "Use Mumps sparse direct solver")
SET(WITH_Hypre FALSE CACHE BOOL "Use Hypre linear algebra library")
SET(WITH_ELMERGUI FALSE CACHE BOOL "Include ElmerGUI")
SET(WITH_ElmerIce FALSE CACHE BOOL "Include ElmerIce")
SET(WITH_ELMERGUITESTER FALSE CACHE BOOL "Include ElmerGUI tester")
SET(WITH_ELMERGUILOGGER FALSE CACHE BOOL "Include ElmerGUI logger")
SET(ELMER_INSTALL_LIB_DIR "lib/elmersolver" CACHE PATH "Location of elmer shared libraries relative to prefix.")
SET(WITH_Trilinos FALSE CACHE BOOL "Use Trilinos")
SET(WITH_FETI4I FALSE CACHE BOOL "Use FETI4I")
SET(WITH_ELMERPOST FALSE CACHE BOOL "Include ElmerPost (DEPRECATED)")
SET(WITH_CONTRIB FALSE CACHE BOOL "Include contributed solvers")

MARK_AS_ADVANCED(WITH_ELMERPOST)

if("${CMAKE_VERSION}" VERSION_GREATER 2.8.12)
  CMAKE_POLICY(SET CMP0022 OLD)
ENDIF()

MARK_AS_ADVANCED(ELMER_INSTALL_LIB_DIR)

ENABLE_TESTING()

SET(ELMER_FEM_MAJOR_VERSION 8)
SET(ELMER_FEM_MINOR_VERSION 3)

SET(ELMER_FEM_VERSION
  ${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION})

IF(NOT(ELMER_FEM_REVISION))
  SET(RESET_ELMER_REVISION TRUE CACHE BOOL "")
  MARK_AS_ADVANCED(RESET_ELMER_REVISION)
ENDIF()

IF(RESET_ELMER_REVISION)
  execute_process(
    COMMAND "git" "--no-pager" "log" "-1" "--pretty=format:%h"
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE ELMER_FEM_REVISION_
    RESULT_VARIABLE RUN_RETURN_VALUE)
  IF(RUN_RETURN_VALUE EQUAL 0)
    SET(ELMER_FEM_REVISION ${ELMER_FEM_REVISION_} CACHE STRING "" FORCE)
    MARK_AS_ADVANCED(ELMER_FEM_REVISION)
  ENDIF()
ENDIF()
    
SET(VERSION ${ELMER_FEM_VERSION})
SET(REVISION ${ELMER_FEM_REVISION})
SET(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
STRING(TIMESTAMP ELMER_FEM_COMPILATIONDATE "%Y-%m-%d")
SET(COMPILATIONDATE ${ELMER_FEM_COMPILATIONDATE})

# Check that the gfortran version is >= 4.8
IF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  INCLUDE(testGFortranVersion)
  IF(NOT CMAKE_Fortran_COMPILER_GNU_VERSION_OK)
    MESSAGE(FATAL_ERROR "GNU Fortran version is too old, should be at least 4.8")
  ENDIF()
ENDIF()

IF(WITH_OpenMP)
  # Advanced properties
  MARK_AS_ADVANCED(
    OpenMP_C_FLAGS
    OpenMP_Fortran_FLAGS
    OpenMP_CXX_FLAGS
    )
  
  FIND_PACKAGE(OpenMP REQUIRED)
  # Add OpenMP flags to compilation flags
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

IF(WITH_MPI)
  # Advanced properties
  MARK_AS_ADVANCED(
    MPI_EXTRA_LIBRARY
    MPI_LIBRARY
    )


  FIND_PACKAGE(MPI REQUIRED)
  # Workaround if MPI is found but MPI_FOUND is not set (windows+msmpi+findmpi)
  SET(MPI_FOUND TRUE CACHE BOOL "")
  MARK_AS_ADVANCED(MPI_FOUND)
  INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
  # MPI libraries may be needed for testing
  SET(CMAKE_C_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
  SET(CMAKE_CXX_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})
  SET(CMAKE_Fortran_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})
ENDIF(WITH_MPI)

IF(WITH_MKL)
  # Find MKL
  FIND_PACKAGE(MKL REQUIRED)
  SET(HAVE_MKL TRUE)
  MARK_AS_ADVANCED(HAVE_MKL)
  # Set BLAS and LAPACK to point to MKL
  SET(BLAS_LIBRARIES ${MKL_BLAS_LIBRARIES} CACHE FILEPATH "")
  SET(LAPACK_LIBRARIES ${MKL_LAPACK_LIBRARIES} CACHE FILEPATH "")
  IF(MPI_FOUND)
    SET(SCALAPACK_LIBRARIES ${MKL_SCALAPACK_LIBRARIES})
    # Versions >11.1 of MKL include Cluster PARDISO direct solver
    IF(WITH_OpenMP AND MKL_CPARDISO_FOUND)
      SET(HAVE_CPARDISO TRUE)
      MARK_AS_ADVANCED(HAVE_CPARDISO)
    ENDIF()
  ENDIF()
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_Fortran_FLAGS}")
ENDIF()

# MARK_AS_ADVANCED(WITH_Trilinos)
IF(WITH_Trilinos)
  FIND_PACKAGE(ML REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/ML ${TRILINOS_PATH}) 
  FIND_PACKAGE(Trilinos REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH})
  FIND_PACKAGE(Belos REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/Belos ${TRILINOS_PATH})
  IF(ML_DIR)
    ADD_DEFINITIONS(-DHAVE_TRILINOS)
  ENDIF()
ENDIF()

# MPI and LAPACK
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
# BLAS and LAPACK libraries may be needed for testing
SET(CMAKE_REQUIRED_LIBRARIES 
  ${CMAKE_REQUIRED_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

IF(MPI_FOUND)
  # Mumps
  IF(WITH_Mumps)
    FIND_PACKAGE(Mumps REQUIRED)
  ENDIF()
  # Hypre
  IF(WITH_Hypre)
    FIND_PACKAGE(Hypre REQUIRED OPTIONAL_COMPONENTS Euclid parcsr_ls ParaSails IJ_mv QUIET)
  ENDIF()

  # Permon
  IF(WITH_FETI4I)
    FIND_PACKAGE(FETI4I REQUIRED)
  ENDIF()
ENDIF()

IF(Mumps_FOUND)
  SET(HAVE_MUMPS TRUE)
  MARK_AS_ADVANCED(HAVE_MUMPS)
  # LINK_DIRECTORIES(${MUMPS_LIBRARIES})
ENDIF()

IF(Hypre_FOUND)
  SET(HAVE_HYPRE TRUE)
  MARK_AS_ADVANCED(HAVE_HYPRE)
  INCLUDE_DIRECTORIES(${Hypre_INCLUDE_DIR})
  # LINK_DIRECTORIES(${Hypre_LIBRARIES})
ENDIF()

IF(FETI4I_FOUND)
  SET(HAVE_FETI4I TRUE)
  MARK_AS_ADVANCED(HAVE_FETI4I)
  ADD_DEFINITIONS(-DHAVE_FETI4I)
ENDIF()

# Check if Fortran compiler supports procedure pointer
INCLUDE(testProcedurePointer)
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_PROCEDUREPOINTER)
  MESSAGE(FATAL_ERROR "Fortran compiler does not seem to support the PROCEDURE statement.")
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
  SET(ELMERLIB_RPATH_STRING "\$ORIGIN/")
ELSEIF(APPLE)                               
  SET(ELMERSOLVER_RPATH_STRING "@loader_path/../${ELMER_INSTALL_LIB_DIR}")
  SET(ELMERLIB_RPATH_STRING "@loader_path/")
ENDIF()

# Uncomment these to use rpath with linked libraries
# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) 
# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# Remove system directories (from the CMake guide)
# LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
# IF("${isSystemDir}" STREQUAL "-1")
#    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# ENDIF("${isSystemDir}" STREQUAL "-1")

ADD_SUBDIRECTORY(matc)

IF(WITH_ElmerIce)
  MESSAGE(STATUS "Adding optional package ElmerIce")
  ADD_SUBDIRECTORY(elmerice)
ENDIF(WITH_ElmerIce)

ADD_SUBDIRECTORY(umfpack)
ADD_SUBDIRECTORY(fhutiter)
ADD_SUBDIRECTORY(meshgen2d)
ADD_SUBDIRECTORY(fem)
ADD_SUBDIRECTORY(mathlibs)
ADD_SUBDIRECTORY(elmergrid)
ADD_SUBDIRECTORY(license_texts)

IF(WITH_ELMERGUI)
    MESSAGE(STATUS "  Building ElmerGUI")
    MESSAGE(STATUS "------------------------------------------------")
    ADD_SUBDIRECTORY(ElmerGUI)
ENDIF(WITH_ELMERGUI)

IF(WITH_ELMERGUITESTER)
  ADD_SUBDIRECTORY(ElmerGUItester)
ENDIF(WITH_ELMERGUITESTER)

IF(WITH_ELMERGUILOGGER)
  ADD_SUBDIRECTORY(ElmerGUIlogger)
ENDIF(WITH_ELMERGUILOGGER)

#INCLUDE(FeatureSummary)
#FEATURE_SUMMARY(WHAT ENABLED_FEATURES
#  INCLUDE_QUIET_PACKAGES
#  DESCRIPTION "Enabled Features:"
#  VAR enabledFeaturesText)
#MESSAGE(STATUS "${enabledFeaturesText}")

# 
IF(CMAKE_BUILD_TYPE)
  STRING(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPCASE)
  MARK_AS_ADVANCED(CMAKE_BUILD_TYPE_UPCASE)
ENDIF(CMAKE_BUILD_TYPE)

MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  BLAS library:   " "${BLAS_LIBRARIES}")
MESSAGE(STATUS "  LAPACK library: " "${LAPACK_LIBRARIES}")
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  Fortran compiler:        " "${CMAKE_Fortran_COMPILER}")
MESSAGE(STATUS "  Fortran flags:           " "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  C compiler:              " "${CMAKE_C_COMPILER}")
MESSAGE(STATUS "  C flags:                 " "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  CXX compiler:            " "${CMAKE_CXX_COMPILER}")
MESSAGE(STATUS "  CXX flags:               " "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
MESSAGE(STATUS "------------------------------------------------")
IF(WITH_MPI)
  MESSAGE(STATUS "  MPI Fortran:             " "${MPI_Fortran_FOUND}")
  MESSAGE(STATUS "  MPI Fortran compiler:    " "${MPI_Fortran_COMPILER}")
  MESSAGE(STATUS "  MPI Fortran flags:       " "${MPI_Fortran_COMPILE_FLAGS}")
  MESSAGE(STATUS "  MPI Fortran include dir: " "${MPI_Fortran_INCLUDE_PATH}")
  MESSAGE(STATUS "  MPI Fortran libraries:   " "${MPI_Fortran_LIBRARIES}")
  MESSAGE(STATUS "  MPI Fortran link flags:  " "${MPI_Fortran_LINK_FLAGS}")
  MESSAGE(STATUS "------------------------------------------------")
  MESSAGE(STATUS "  MPI C:             " "${MPI_C_FOUND}")
  MESSAGE(STATUS "  MPI C compiler:    " "${MPI_C_COMPILER}")
  MESSAGE(STATUS "  MPI C flags:       " "${MPI_C_COMPILE_FLAGS}")
  MESSAGE(STATUS "  MPI C include dir: " "${MPI_C_INCLUDE_PATH}")
  MESSAGE(STATUS "  MPI C libraries:   " "${MPI_C_LIBRARIES}")
  MESSAGE(STATUS "  MPI C flags:       " "${MPI_C_LINK_FLAGS}")
  MESSAGE(STATUS "------------------------------------------------")
  IF(WITH_Mumps)
    MESSAGE(STATUS "  Mumps:             " "${Mumps_FOUND}")
    MESSAGE(STATUS "  Mumps include:     " "${Mumps_INCLUDE_DIR}")
    MESSAGE(STATUS "  Mumps libraries:   " "${Mumps_LIBRARIES}")
    MESSAGE(STATUS "------------------------------------------------")
  ENDIF(WITH_Mumps)
  IF(WITH_Hypre)
    MESSAGE(STATUS "  Hypre:             " "${Hypre_FOUND}")
    MESSAGE(STATUS "  Hypre include:     " "${Hypre_INCLUDE_DIR}")
    MESSAGE(STATUS "  Hypre libraries:   " "${Hypre_LIBRARIES}")
    MESSAGE(STATUS "------------------------------------------------")
  ENDIF(WITH_Hypre)
ENDIF(WITH_MPI)

IF(WITH_Trilinos)
  MESSAGE(STATUS "  Trilinos:           " "${Trilinos_FOUND}")
  MESSAGE(STATUS "  Trilinos_DIR:       " "${Trilinos_DIR}")
  MESSAGE(STATUS "  Trilinos_LIBRARIES: " "${Trilinos_LIBRARIES}")
  MESSAGE(STATUS "------------------------------------------------")
ENDIF(WITH_Trilinos)

IF(WITH_CONTRIB)
  MESSAGE(STATUS "  Building contributed solvers")
  MESSAGE(STATUS "------------------------------------------------")
ENDIF(WITH_CONTRIB)

IF(WITH_ELMERGUITESTER)
  MESSAGE(STATUS "  Building ElmerGUI tester")
  MESSAGE(STATUS "------------------------------------------------")
ENDIF(WITH_ELMERGUITESTER)

IF(WITH_ELMERGUILOGGER)
  MESSAGE(STATUS "  Building ElmerGUI logger")
  MESSAGE(STATUS "------------------------------------------------")
ENDIF(WITH_ELMERGUILOGGER)

IF(WITH_ELMERPOST)
  MESSAGE(WARNING "  Building ElmerPost")
  MESSAGE(STATUS "------------------------------------------------")
  ADD_SUBDIRECTORY(post)
ENDIF(WITH_ELMERPOST)

# Packaging
if(NOT BYPASS_CPACK)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cpack/ElmerCPack.cmake)
endif(NOT BYPASS_CPACK)
