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

if(APPLE)
  # message("you need to have gcc-gfrotran installed using HomeBrew")
  # set(CMAKE_C_COMPILER "/usr/bin/gcc")
  # set(CMAKE_CXX_COMPILER "/usr/bin/g++")
  set(CMAKE_C_COMPILER "/usr/local/bin/gcc-10")
  set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-10")
  set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran")
  # set(BLA_VENDOR "OpenBLAS")
  # option(HUNTER_ENABLED "Enable Hunter package manager support" OFF)
  # set (CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)
  # set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
endif()





IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type (Release, Debug, RelWithDebugInfo, MinSizeRel)")
ENDIF()

PROJECT(Elmer Fortran C CXX)

# CMAKE_VERSION seems to require this in minimum
CMAKE_MINIMUM_REQUIRED(VERSION 1.0.2)

IF(APPLE)
  SET(CMAKE_MACOSX_RPATH 1) 
    # on OSX FindOpenGL uses framework version of OpenGL, but we need X11 version
    FIND_PATH(GLX_INCLUDE_DIR GL/glx.h
      PATHS /opt/X11/include /usr/X11/include /usr/X11R6/include)
    find_library(X11_GL_LIB GL
      PATH /opt/X11/lib /usr/X11/lib /usr/X11R6/lib)
# FIND_LIBRARY(OPENGL_gl_LIBRARY GL
#        PATHS /opt/X11/lib /usr/X11/lib /usr/X11R6/lib)
# FIND_LIBRARY(OPENGL_glu_LIBRARY GLU
#        PATHS /opt/X11/lib /usr/X11/lib /usr/X11R6/lib)
# LIST(APPEND LIBS ${OPENGL_gl_LIBRARY})
INCLUDE_DIRECTORIES(${GLX_INCLUDE_DIR})
ENDIF()

#shamelessly borrowed from FreeCAD project: https://github.com/FreeCAD/FreeCAD/blob/master/cMake/FreeCAD_Helpers/SetupPython.cmake
# For building on OS X
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  # If the user doesn't tell us which package manager they're using
  if(NOT DEFINED MACPORTS_PREFIX AND NOT DEFINED HOMEBREW_PREFIX)

    # Try to find MacPorts path
    find_program(MACPORTS_EXECUTABLE port)
    if(EXISTS ${MACPORTS_EXECUTABLE})
      string(REPLACE "/bin/port" ""
      MACPORTS_PREFIX ${MACPORTS_EXECUTABLE})
      message(STATUS "Detected MacPorts install at ${MACPORTS_PREFIX}")
    endif(EXISTS ${MACPORTS_EXECUTABLE})

    # Try to find Homebrew path
    find_program(HOMEBREW_EXECUTABLE brew)
    if(EXISTS ${HOMEBREW_EXECUTABLE})
      string(REPLACE "/bin/brew" ""
      HOMEBREW_PREFIX ${HOMEBREW_EXECUTABLE})
      message(STATUS "Detected Homebrew install at ${HOMEBREW_PREFIX}")
    endif()

  endif(NOT DEFINED MACPORTS_PREFIX AND NOT DEFINED HOMEBREW_PREFIX)

  # In case someone tries to shoot themselves in the foot
  if(DEFINED MACPORTS_PREFIX AND DEFINED HOMEBREW_PREFIX)
    message(SEND_ERROR "Multiple package management systems detected - ")
    message(SEND_ERROR "define either MACPORTS_PREFIX or HOMEBREW_PREFIX")

  # No package manager
  elseif(NOT DEFINED MACPORTS_PREFIX AND NOT DEFINED HOMEBREW_PREFIX)
    message(SEND_ERROR "No package manager detected - install MacPorts or Homebrew")
  endif()
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


# SET(CMAKE_MODULE_PATH "${CMAKE_ROOT}/Modules;${PROJECT_SOURCE_DIR}/cmake/Modules")
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules;${CMAKE_ROOT}/Modules")

SET(CMAKE_Fortran_MODULE_DIRECTORY
  ${PROJECT_BINARY_DIR}/fmodules CACHE PATH "Directory for Fortran modules")

# Defaults to cache variables
SET(WITH_MPI TRUE CACHE BOOL "Use MPI parallelization")
SET(MPI_TEST_MAXPROC 8 CACHE STRING "Maximum number of tasks used in parallel tests")
SET(MPI_TEST_MINPROC 1 CACHE STRING "Minimum number of tasks used in parallel tests")
SET(WITH_OpenMP FALSE CACHE BOOL "Use OpenMP thread and SIMD (if available) parallelization")
SET(WITH_MKL FALSE CACHE BOOL "Use Intel Math Kernel library")
SET(WITH_Mumps FALSE CACHE BOOL "Use Mumps sparse direct solver")
SET(WITH_Hypre FALSE CACHE BOOL "Use Hypre linear algebra library")
SET(WITH_ELMERGUI FALSE CACHE BOOL "Include ElmerGUI")
SET(WITH_ElmerIce FALSE CACHE BOOL "Include ElmerIce")
SET(WITH_XIOS FALSE CACHE BOOL "Include XIOS")
SET(WITH_ELMERGUITESTER FALSE CACHE BOOL "Include ElmerGUI tester")
SET(WITH_ELMERGUILOGGER FALSE CACHE BOOL "Include ElmerGUI logger")
SET(ELMER_INSTALL_LIB_DIR "lib/elmersolver" CACHE PATH "Location of elmer shared libraries relative to prefix.")
SET(WITH_Trilinos FALSE CACHE BOOL "Use Trilinos")
SET(WITH_FETI4I FALSE CACHE BOOL "Use FETI4I")
SET(WITH_ELMERPOST FALSE CACHE BOOL "Include ElmerPost (DEPRECATED)")
SET(WITH_CONTRIB FALSE CACHE BOOL "Include contributed solvers")

SET(WITH_LUA FALSE CACHE BOOL "Include LUA extensions.")
SET(WITH_MMG FALSE CACHE BOOL "Compile with MMG library.")

SET(WITH_Zoltan FALSE CACHE BOOL "Link in Zoltan mesh repartitioning library.")

# Let user select if the contiguous attribute is used or not
SET(USE_CONTIGUOUS FALSE CACHE BOOL "Use contiguous attribute")

MARK_AS_ADVANCED(WITH_ELMERPOST)

# if("${CMAKE_VERSION}" VERSION_GREATER 2.8.12)
#   CMAKE_POLICY(SET CMP0022 OLD)
# ENDIF()

OPTION (USE_NEW_POLICY_CMP0057 "Option to enable new policy CMP0057" ON)
IF (USE_NEW_POLICY_CMP0057)
  CMAKE_POLICY(SET CMP0057 NEW)
ENDIF()

MARK_AS_ADVANCED(ELMER_INSTALL_LIB_DIR)

ENABLE_TESTING()

SET(ELMER_FEM_MAJOR_VERSION 9)
SET(ELMER_FEM_MINOR_VERSION 0)

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

# Check that the gfortran version is >= 7
IF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  INCLUDE(testGFortranVersion)
  IF(NOT CMAKE_Fortran_COMPILER_GNU_VERSION_OK)
    MESSAGE(FATAL_ERROR "GNU Fortran major version is too old, should be at least 7")
  ENDIF()
ENDIF()

IF(WITH_LUA)
  IF(NOT(${USE_SYSTEM_LUA}))
    ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/contrib/lua-5.1.5/)
    set(LUA_LIBRARIES luaelmer CACHE STRING "")
    set(LUA_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/contrib/lua-5.1.5/src CACHE STRING "")
  ELSE()
    FIND_PACKAGE(Lua51 REQUIRED)
  ENDIF()
  SET(HAVE_LUA TRUE)
  MARK_AS_ADVANCED(HAVE_LUA)
ENDIF()

IF(WITH_EXTOPTIM)
  SET(HAVE_EXTOPTIM TRUE)
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
  # if(APPLE)
  #   if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  #     set(OpenMP_C "${CMAKE_C_COMPILER}")
  #     set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
  #     set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
  #     set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
  #     set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
  #     set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
  #   endif()
  #   if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  #     set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
  #     set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
  #     set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  #     set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  #     set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  #     set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  #   endif()
  # else()
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  # endif()

 
  
  # Test compiler support for OpenMP 4.0 features used
  INCLUDE(testOpenMP40)
  IF(CMAKE_Fortran_COMPILER_SUPPORTS_OPENMP40)
    SET(HAVE_OPENMP40 TRUE)
    MARK_AS_ADVANCED(HAVE_OPENMP40)
    # Test compiler support for advanced OpenMP 4.5 features used
    INCLUDE(testOpenMP45)
    IF(CMAKE_Fortran_COMPILER_SUPPORTS_OPENMP45)
      SET(HAVE_OPENMP45 TRUE)
      MARK_AS_ADVANCED(HAVE_OPENMP45)
    ENDIF()
  ENDIF()
ENDIF()

# Get rid of the annoying rank mismatch warning
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  IF(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.9)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
#    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")
  ENDIF()
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

# if(APPLE)
#   # message(STATUS "Trying to find Accelerate library")
#   # find_library(BLAS_LIBRARIES Accelerate HINTS /System/Library/Frameworks/Accelerate.framework)
#   # MARK_AS_ADVANCED(BLAS_LIBRARIES)
#   # set(BLAS_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -framework accelerate")
#   # find_path(LAPACK_INCLUDE_DIRS
#   #   NAMES clapack.h
#   #   PATHS ${BLAS_INCLUDE_DIRS}
#   #   NO_DEFAULT_PATH
#   # )
#   # set(LAPACK_HEADER ${LAPACK_INCLUDE_DIRS}/clapack.h)
#   # set(LAPACK_SUFFIX "_")
#   # set(LAPACK_PREFIX)

#   # message(STATUS "Using Accelerate framework on macOS")
#   message("using the local blas library")
#   # option(BLAS_VERBOSE "Print some additional information during BLAS libraries detection" ON)
#   set(BLAS_DIR "./mathlibs/src/blas/")
# set(BLA_VENDOR "Generic")

  # set(ACCELERATE_LIBS "-framework Accelerate")
  # BLAS_Accelerate_LIBRARY:FILEPATH=/System/Library/Frameworks/Accelerate.framework


# else()
  IF(WITH_MKL)
    # Find MKL
    FIND_PACKAGE(MKL REQUIRED)
    SET(HAVE_MKL TRUE)
    MARK_AS_ADVANCED(HAVE_MKL)
    # Set BLAS and LAPACK to point to MKL
    SET(BLAS_LIBRARIES ${MKL_BLAS_LIBRARIES} CACHE FILEPATH "")
    SET(LAPACK_LIBRARIES ${MKL_LAPACK_LIBRARIES} CACHE FILEPATH "")
    IF(MKL_SCALAPACK_LIBRARIES_FOUND)
      SET(SCALAPACK_LIBRARIES ${MKL_SCALAPACK_LIBRARIES})
      # Versions >11.1 of MKL include Cluster PARDISO direct solver
      IF(MKL_CPARDISO_FOUND)
        SET(HAVE_CPARDISO TRUE)
        MARK_AS_ADVANCED(HAVE_CPARDISO)
      ENDIF()
    ENDIF()
    INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_Fortran_FLAGS}")
  ELSE()
    IF (NOT BLAS_LIBRARIES OR NOT LAPACK_LIBRARIES)
      # If BLAS/LAPACK not explicitly set, attempt to find MKL BLAS/LAPACK
      FIND_PACKAGE(MKL QUIET)
      IF (MKL_BLAS_LIBRARIES_FOUND AND MKL_LAPACK_LIBRARIES_FOUND)
        # Set BLAS and LAPACK to point to MKL, in order to work around
        # CMake FindBLAS finding wrong core library for GNU Fortran
        SET(BLAS_LIBRARIES ${MKL_BLAS_LIBRARIES} CACHE FILEPATH "")
        SET(LAPACK_LIBRARIES ${MKL_LAPACK_LIBRARIES} CACHE FILEPATH "")
        INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_Fortran_FLAGS}")
      ENDIF()
    ENDIF()
  ENDIF()
# endif(APPLE)

# MARK_AS_ADVANCED(WITH_Trilinos)
IF(WITH_Trilinos)
  FIND_PACKAGE(ML REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/ML ${TRILINOS_PATH})
  FIND_PACKAGE(Trilinos REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH})
  FIND_PACKAGE(Belos REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/Belos ${TRILINOS_PATH})
  IF(ML_DIR)
    ADD_DEFINITIONS(-DHAVE_TRILINOS)
  ENDIF()
ENDIF()

IF(WITH_Zoltan)
  # ZOLTAN Remeshing
  MESSAGE(STATUS "------------------------------------------------")

  SET(USE_SYSTEM_ZOLTAN FALSE CACHE BOOL "Use system Zoltan library")
  IF(NOT ${USE_SYSTEM_ZOLTAN})
    MESSAGE(STATUS "Building Zoltan mesh repartitioner")
    ADD_SUBDIRECTORY(contrib/Zoltan_v3.83)
    SET(HAVE_ZOLTAN TRUE CACHE BOOL "Has Zoltan tools for mesh (re)partitioning" )
    MARK_AS_ADVANCED(HAVE_ZOLTAN)
    MESSAGE(STATUS "  Zoltan_LIB:      " "${ZOLTAN_LIBRARY}")
    MESSAGE(STATUS "------------------------------------------------")
  ELSE()
    FIND_PACKAGE(Zoltan)

    IF(ZOLTAN_FOUND)

      SET(HAVE_ZOLTAN TRUE CACHE BOOL "Has Zoltan tools for mesh (re)partitioning" )
      GET_FILENAME_COMPONENT(ZOLTAN_LIBDIR ${ZOLTAN_LIBRARY} DIRECTORY)

      MARK_AS_ADVANCED(HAVE_ZOLTAN)
      MESSAGE(STATUS "  Zoltan:           " "${ZOLTAN_FOUND}")
      MESSAGE(STATUS "  Zoltan_INC:       " "${ZOLTAN_INCLUDE_DIR}")
      MESSAGE(STATUS "  Zoltan_LIB:      " "${ZOLTAN_LIBRARY}")
      MESSAGE(STATUS "  Zoltan_LIBDIR:      " "${ZOLTAN_LIBDIR}")

      INCLUDE_DIRECTORIES(${ZOLTAN_INCLUDE_DIR})
      MESSAGE(STATUS "------------------------------------------------")

    ELSE()
      MESSAGE(STATUS "  Library not found: >ZOLTAN_FOUND< ")
      MESSAGE(STATUS "    Missing: >ZOLTAN_INCLUDE_DIR< , >ZOLTAN_LIBRARY< for mesh (re)partitioning")
    ENDIF(ZOLTAN_FOUND)

    MESSAGE(STATUS "------------------------------------------------")
  ENDIF()
ENDIF(WITH_Zoltan)


IF(WITH_NETCDF)
  # ---------------------- #
  # -- NETCDF libraries -- #
  #----------------------- #
  MESSAGE(STATUS "------------------------------------------------")
  MESSAGE(STATUS "Elmer/Ice package: Looking for [NETCDF] & [NETCDF Fortran] libraries")
  
  
  FIND_PACKAGE(NETCDF)
  
  IF(NETCDF_FOUND)
  
    SET(HAVE_NETCDF TRUE CACHE BOOL "Has NETCDF for Elmer/Ice Output" )
    GET_FILENAME_COMPONENT(NETCDF_LIBDIR ${NETCDF_LIBRARY} DIRECTORY)
    GET_FILENAME_COMPONENT(NETCDFF_LIBDIR ${NETCDFF_LIBRARY} DIRECTORY)
  
    MARK_AS_ADVANCED(HAVE_NETCDF)
    INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIR})
    ADD_DEFINITIONS(-DHAVE_NETCDF)
  
    MESSAGE(STATUS "  NETCDF:           " "${NETCDF_FOUND}")
    MESSAGE(STATUS "  NETCDF_INC:       " "${NETCDF_INCLUDE_DIR}")
    MESSAGE(STATUS "  NETCDF_LIBS:     " "${NETCDF_LIBRARIES}")
  ELSE()
    MESSAGE(STATUS "Library not found: >NETCDF_FOUND< ")
    MESSAGE(WARNING "  \n Missing: >NETCDF_INCLUDE_DIR< , >NETCDF_LIBRARY<,  >NETCDFF_LIBRARY< \n some functionalities will be disabled")
  ENDIF(NETCDF_FOUND)
  
  IF(NETCDF_FOUND)
    # ---------------------- #
    # -- HDF5 libraries -- #
    #----------------------- #
    MESSAGE(STATUS "------------------------------------------------")
    MESSAGE(STATUS "Elmer/Ice package: Looking for [HDF5] & [HDF5 Fortran] libraries")
  
    FIND_PACKAGE(PHDF5)
    
    IF(PHDF5_FOUND)
    
      SET(HAVE_HDF5 TRUE CACHE BOOL "Has HDF5 for Elmer/Ice Output" )
      GET_FILENAME_COMPONENT(PHDF5_LIBDIR ${PHDF5_LIBRARY} DIRECTORY)
      GET_FILENAME_COMPONENT(PHDF5HL_LIBDIR ${PHDF5HL_LIBRARY} DIRECTORY)
    
      MARK_AS_ADVANCED(HAVE_HDF5)
      INCLUDE_DIRECTORIES(${PHDF5_INCLUDE_DIR})
      ADD_DEFINITIONS(-DHAVE_HDF5)
    
      MESSAGE(STATUS "  HDF5:           " "${PHDF5_FOUND}")
      MESSAGE(STATUS "  HDF5_INC:       " "${PHDF5_INCLUDE_DIR}")
      MESSAGE(STATUS "  HDF5_LIBS:     " "${PHDF5_LIBRARIES}")
    ELSE()
      MESSAGE(STATUS "Library not found: >HDF5_FOUND< ")
      MESSAGE(WARNING "  \n Missing: >HDF5_INCLUDE_DIR< , >HDF5_LIBRARY<,  >HDF5F_LIBRARY< \n some functionalities will be disabled")
    ENDIF(PHDF5_FOUND)
  ENDIF(NETCDF_FOUND)
  
ENDIF(WITH_NETCDF)

# MMG Remeshing
IF(WITH_MMG)
  MESSAGE(STATUS "------------------------------------------------")
  MESSAGE(STATUS "Mesh adaptation 2D/3D looking for [Mmg] tools ")
  FIND_PACKAGE(MMG)
  # INCLUDE("cmake/Modules/FindMMG.cmake")
  # INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR})
  # TARGET_LINK_LIBRARIES( ${YOUR_TARGET} ${MMG_LIBRARY})
  IF(MMG_FOUND)

    SET(HAVE_MMG TRUE CACHE BOOL "Has Mmg tools for Elmer/Ice Mesh adaptation" )
    GET_FILENAME_COMPONENT(MMG_LIBDIR ${MMG_LIBRARY} DIRECTORY)

    MARK_AS_ADVANCED(HAVE_MMG)
    MESSAGE(STATUS "  Mmg:           " "${MMG_FOUND}")
    MESSAGE(STATUS "  Mmg_INC:       " "${MMG_INCLUDE_DIR}")
    MESSAGE(STATUS "  Mmg_LIB:      " "${MMG_LIBRARY}")
    MESSAGE(STATUS "  Mmg_LIBDIR:      " "${MMG_LIBDIR}")
    MESSAGE(STATUS "Compile MMG2DSolver/MMG3DSolver")
    ADD_DEFINITIONS(-DHAVE_MMG)

    INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR})
    MESSAGE(STATUS "------------------------------------------------")

  ELSE()
    MESSAGE(STATUS "  Library not found: >MMG_FOUND< ")
    MESSAGE(STATUS "    Missing: >MMG_INCLUDE_DIR< , >MMG_LIBRARY<, to compile MMG3DSolver")
  ENDIF(MMG_FOUND)

  MESSAGE(STATUS "------------------------------------------------")

  IF (WITH_PARMMG) 
    FIND_PACKAGE(PARMMG)
    IF(PARMMG_FOUND)

      SET(HAVE_PARMMG TRUE CACHE BOOL "Has ParMmg tools for Elmer/Ice Mesh adaptation" )
      GET_FILENAME_COMPONENT(PARMMG_LIBDIR ${PARMMG_LIBRARY} DIRECTORY)
      
      MARK_AS_ADVANCED(HAVE_PARMMG)
      MESSAGE(STATUS "  ParMmg:           " "${PARMMG_FOUND}")
      MESSAGE(STATUS "  ParMmg_INC:       " "${PARMMG_INCLUDE_DIR}")
      MESSAGE(STATUS "  ParMmg_LIB:      " "${PARMMG_LIBRARY}")
      MESSAGE(STATUS "  ParMmg_LIBDIR:      " "${PARMMG_LIBDIR}")
      #MESSAGE(STATUS "Compile MMG2DSolver/MMG3DSolver")
      ADD_DEFINITIONS(-DHAVE_PARMMG)

      INCLUDE_DIRECTORIES(${PARMMG_INCLUDE_DIR})
      MESSAGE(STATUS "------------------------------------------------")
    ELSE()
      MESSAGE(STATUS "  Library not found: >PARMMG_FOUND< ")
      MESSAGE(STATUS "    Missing: >PARMMG_INCLUDE_DIR< , >PARMMG_LIBRARY<, to compile PARMMGSolver")
    ENDIF(PARMMG_FOUND)
  ELSE()
    SET(PARMMG_FOUND FALSE)
    MESSAGE(STATUS "  Skipping ParMMG")
  # INCLUDE("cmake/Modules/FindParMMG.cmake")
  # INCLUDE_DIRECTORIES(${PARMMG_INCLUDE_DIR})
  # TARGET_LINK_LIBRARIES( ${YOUR_TARGET} ${PARMMG_LIBRARY})
  #  SET(HAVE_PARMMG FALSE)
  ENDIF(WITH_PARMMG) 
ELSE()
  SET(HAVE_MMG FALSE)
  SET(HAVE_PARMMG FALSE)
ENDIF()

# XIOS
IF (WITH_XIOS)
  MESSAGE(STATUS "------------------------------------------------")
  MESSAGE(STATUS "  Looking for XIOS ")

  FIND_PACKAGE(XIOS)

  IF(XIOS_FOUND)
    SET(HAVE_XIOS TRUE)
    MARK_AS_ADVANCED(HAVE_XIOS)
    INCLUDE_DIRECTORIES(${XIOS_INCLUDE_DIR})
    LINK_DIRECTORIES(${XIOS_LIBRARIES})
    ADD_DEFINITIONS(-DHAVE_XIOS)
    MESSAGE(STATUS "  XIOS:           " "${XIOS_FOUND}")
    MESSAGE(STATUS "  XIOS_INC:       " "${XIOS_INCLUDE_DIR}")
    MESSAGE(STATUS "  XIOS_LIB :      " "${XIOS_LIBRARIES}")
  ELSE()
    MESSAGE(STATUS "  XIOS not found ")
  ENDIF()

  MESSAGE(STATUS "------------------------------------------------")
ENDIF()

MESSAGE(STATUS "------------------------------------------------")

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
    FIND_PACKAGE(Hypre REQUIRED OPTIONAL_COMPONENTS Euclid parcsr_ls ParaSails IJ_mv FEI core QUIET)
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
IF(USE_CONTIGUOUS AND CMAKE_Fortran_COMPILER_SUPPORTS_CONTIGUOUS)
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

IF(HAVE_LUA)
  MESSAGE(STATUS "  Building Lua")
  MESSAGE(STATUS "------------------------------------------------")
ENDIF(HAVE_LUA)

IF(HAVE_EXTOPTIM)
  MESSAGE(STATUS "  Building with external optimization routines")
  MESSAGE(STATUS "------------------------------------------------")
ENDIF(HAVE_EXTOPTIM)



# Packaging
if(NOT BYPASS_CPACK)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cpack/ElmerCPack.cmake)
endif(NOT BYPASS_CPACK)
