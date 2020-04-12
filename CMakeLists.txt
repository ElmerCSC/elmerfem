#
# CMake configuration script for Elmer
#
# Authors: Sami Ilvonen, CSC - IT Center for Science, Ltd. Phil Weir, NUMA
# Engineering Services, Ltd. Mikko Byckling, CSC - IT Center for Science Ltd.
# Juhani Kataja, CSC - IT Center for Science Ltd.
#
# First public beta release 11th of October, 2013
#

if(APPLE)
  # message("you need to have gcc-gfrotran installed using HomeBrew")
  set(CMAKE_C_COMPILER "/usr/local/bin/gcc-9")
  set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-9")
  set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran")
  # set(BLA_VENDOR "OpenBLAS") option(HUNTER_ENABLED "Enable Hunter package
  # manager support" OFF) set (CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL
  # "" FORCE) set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
endif()

set(CMAKE_LEGACY_CYGWIN_WIN32 0)

# if(MSYS) set(CMAKE_C_COMPILER "/mingw64/bin/gcc") set(CMAKE_CXX_COMPILER
# "/mingw64/bin/g++") set(CMAKE_Fortran_COMPILER "/mingw64/bin/gfortran")
# set (CMAKE_GENERATOR "MSYS Makefiles" CACHE INTERNAL "" FORCE) #
# set(CMAKE_MAKE_PROGRAM /mingw64/bin/mingw32-make CACHE INTERNAL "") #
# set(CMAKE_MAKE_PROGRAM /mingw64/bin/mingw32-make CACHE FILEPATH "" FORCE)
# endif(MSYS)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE STRING "Build type (Release, Debug, RelWithDebugInfo, MinSizeRel)")
endif()

project(Elmer Fortran C CXX)
cmake_minimum_required(VERSION 2.8.9)

if(APPLE)
  set(CMAKE_MACOSX_RPATH 1)
  # on OSX FindOpenGL uses framework version of OpenGL, but we need X11 version
  find_path(GLX_INCLUDE_DIR GL/glx.h PATHS /opt/X11/include /usr/X11/include
                                           /usr/X11R6/include)
  find_library(X11_GL_LIB GL PATH /opt/X11/lib /usr/X11/lib /usr/X11R6/lib)
  # FIND_LIBRARY(OPENGL_gl_LIBRARY GL PATHS /opt/X11/lib /usr/X11/lib
  # /usr/X11R6/lib) FIND_LIBRARY(OPENGL_glu_LIBRARY GLU PATHS /opt/X11/lib
  # /usr/X11/lib /usr/X11R6/lib) LIST(APPEND LIBS ${OPENGL_gl_LIBRARY})
  include_directories(${GLX_INCLUDE_DIR})
endif()

# shamelessly borrowed from FreeCAD project: https://github.com/FreeCAD/FreeCAD/
# blob/master/cMake/FreeCAD_Helpers/SetupPython.cmake For building on OS X
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  # If the user doesn't tell us which package manager they're using
  if(NOT DEFINED MACPORTS_PREFIX AND NOT DEFINED HOMEBREW_PREFIX)

    # Try to find MacPorts path
    find_program(MACPORTS_EXECUTABLE port)
    if(EXISTS ${MACPORTS_EXECUTABLE})
      string(REPLACE "/bin/port" "" MACPORTS_PREFIX ${MACPORTS_EXECUTABLE})
      message(STATUS "Detected MacPorts install at ${MACPORTS_PREFIX}")
    endif(EXISTS ${MACPORTS_EXECUTABLE})

    # Try to find Homebrew path
    find_program(HOMEBREW_EXECUTABLE brew)
    if(EXISTS ${HOMEBREW_EXECUTABLE})
      string(REPLACE "/bin/brew" "" HOMEBREW_PREFIX ${HOMEBREW_EXECUTABLE})
      message(STATUS "Detected Homebrew install at ${HOMEBREW_PREFIX}")
    endif()

  endif(NOT DEFINED MACPORTS_PREFIX AND NOT DEFINED HOMEBREW_PREFIX)

  # In case someone tries to shoot themselves in the foot
  if(DEFINED MACPORTS_PREFIX AND DEFINED HOMEBREW_PREFIX)
    message(SEND_ERROR "Multiple package management systems detected - ")
    message(SEND_ERROR "define either MACPORTS_PREFIX or HOMEBREW_PREFIX")

    # No package manager
  elseif(NOT DEFINED MACPORTS_PREFIX AND NOT DEFINED HOMEBREW_PREFIX)
    message(
      SEND_ERROR "No package manager detected - install MacPorts or Homebrew")
  endif()
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules"
                      ${CMAKE_ROOT}/Modules)

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/fmodules
    CACHE PATH "Directory for Fortran modules")

# Defaults to cache variables
set(WITH_MPI
    TRUE
    CACHE BOOL "Use MPI parallelization")
set(WITH_OpenMP
    TRUE
    CACHE BOOL "Use OpenMP thread and SIMD (if available) parallelization")
set(WITH_MKL
    FALSE
    CACHE BOOL "Use Intel Math Kernel library")
set(WITH_Mumps
    FALSE
    CACHE BOOL "Use Mumps sparse direct solver")
set(WITH_Hypre
    FALSE
    CACHE BOOL "Use Hypre linear algebra library")
set(WITH_ELMERGUI
    FALSE
    CACHE BOOL "Include ElmerGUI")
set(WITH_ElmerIce
    FALSE
    CACHE BOOL "Include ElmerIce")
set(WITH_ELMERGUITESTER
    FALSE
    CACHE BOOL "Include ElmerGUI tester")
set(WITH_ELMERGUILOGGER
    FALSE
    CACHE BOOL "Include ElmerGUI logger")
set(ELMER_INSTALL_LIB_DIR
    "lib/elmersolver"
    CACHE PATH "Location of elmer shared libraries relative to prefix.")
set(WITH_Trilinos
    FALSE
    CACHE BOOL "Use Trilinos")
set(WITH_FETI4I
    FALSE
    CACHE BOOL "Use FETI4I")
set(WITH_ELMERPOST
    FALSE
    CACHE BOOL "Include ElmerPost (DEPRECATED)")
set(WITH_CONTRIB
    FALSE
    CACHE BOOL "Include contributed solvers")
set(WITH_LUA
    FALSE
    CACHE BOOL "Include LUA extensions.")
set(WITH_Zoltan
    FALSE
    CACHE BOOL "Link in Zoltan mesh repartitioning library.")
set(WITH_MMG
    FALSE
    CACHE BOOL "Include mmg - Surface and volume remesher library.")

# Let user select if the contiguous attribute is used or not
set(USE_CONTIGUOUS
    FALSE
    CACHE BOOL "Use contiguous attribute")

# if("${CMAKE_VERSION}" VERSION_GREATER 2.8.12) cmake_policy(SET CMP0022 OLD)
# endif()

cmake_policy(SET CMP0075 NEW)
cmake_policy(SET CMP0054 NEW)

enable_testing()

set(ELMER_FEM_MAJOR_VERSION 8)
set(ELMER_FEM_MINOR_VERSION 4)
set(ELMER_FEM_VERSION ${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION})

if(NOT (ELMER_FEM_REVISION))
  set(RESET_ELMER_REVISION
      TRUE
      CACHE BOOL "")
  mark_as_advanced(RESET_ELMER_REVISION)
endif()

if(RESET_ELMER_REVISION)
  execute_process(
    COMMAND "git" "--no-pager" "log" "-1" "--pretty=format:%h"
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE ELMER_FEM_REVISION_
    RESULT_VARIABLE RUN_RETURN_VALUE)
  if(RUN_RETURN_VALUE EQUAL 0)
    set(ELMER_FEM_REVISION
        ${ELMER_FEM_REVISION_}
        CACHE STRING "" FORCE)
    mark_as_advanced(ELMER_FEM_REVISION)
  endif()
endif()

set(VERSION ${ELMER_FEM_VERSION})
set(REVISION ${ELMER_FEM_REVISION})
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
string(TIMESTAMP ELMER_FEM_COMPILATIONDATE "%Y-%m-%d")
set(COMPILATIONDATE ${ELMER_FEM_COMPILATIONDATE})

# Check that the gfortran version is >= 4.8
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  include(testGFortranVersion)
  if(NOT CMAKE_Fortran_COMPILER_GNU_VERSION_OK)
    message(
      FATAL_ERROR "GNU Fortran version is too old, should be at least 4.8")
  endif()
endif()

if(WITH_LUA)
  if(NOT (${USE_SYSTEM_LUA}))
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/contrib/lua-5.1.5/)
    set(LUA_LIBRARIES
        luaelmer
        CACHE STRING "")
    set(LUA_INCLUDE_DIR
        ${PROJECT_SOURCE_DIR}/contrib/lua-5.1.5/src
        CACHE STRING "")
  else()
    find_package(Lua51 REQUIRED)
  endif()
  set(HAVE_LUA TRUE)
  mark_as_advanced(HAVE_LUA)
endif()

if(WITH_OpenMP)
  # Advanced properties
  mark_as_advanced(OpenMP_C_FLAGS OpenMP_Fortran_FLAGS OpenMP_CXX_FLAGS)

  find_package(OpenMP REQUIRED)
  # Add OpenMP flags to compilation flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

  # Test compiler support for OpenMP 4.0 features used
  include(testOpenMP40)
  if(CMAKE_Fortran_COMPILER_SUPPORTS_OPENMP40)
    set(HAVE_OPENMP40 TRUE)
    mark_as_advanced(HAVE_OPENMP40)
    # Test compiler support for advanced OpenMP 4.5 features used
    include(testOpenMP45)
    if(CMAKE_Fortran_COMPILER_SUPPORTS_OPENMP45)
      set(HAVE_OPENMP45 TRUE)
      mark_as_advanced(HAVE_OPENMP45)
    endif()
  endif()
endif()

if(WITH_MPI)
  set(MPI_TEST_MAXPROC
      8
      CACHE STRING "Maximum number of tasks used in parallel tests")
  set(MPI_TEST_MINPROC
      1
      CACHE STRING "Minimum number of tasks used in parallel tests")

  # Advanced properties
  mark_as_advanced(MPI_EXTRA_LIBRARY MPI_LIBRARY)

  find_package(MPI REQUIRED)
  # Workaround if MPI is found but MPI_FOUND is not set (windows+msmpi+findmpi)
  set(MPI_FOUND
      TRUE
      CACHE BOOL "")
  mark_as_advanced(MPI_FOUND)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  # MPI libraries may be needed for testing
  set(CMAKE_C_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
  set(CMAKE_CXX_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})
  set(CMAKE_Fortran_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})

  message(STATUS "  MPI Fortran:             " "${MPI_Fortran_FOUND}")
  message(STATUS "  MPI Fortran compiler:    " "${MPI_Fortran_COMPILER}")
  message(STATUS "  MPI Fortran flags:       " "${MPI_Fortran_COMPILE_FLAGS}")
  message(STATUS "  MPI Fortran include dir: " "${MPI_Fortran_INCLUDE_PATH}")
  message(STATUS "  MPI Fortran libraries:   " "${MPI_Fortran_LIBRARIES}")
  message(STATUS "  MPI Fortran link flags:  " "${MPI_Fortran_LINK_FLAGS}")
  message(STATUS "------------------------------------------------")
  message(STATUS "  MPI C:             " "${MPI_C_FOUND}")
  message(STATUS "  MPI C compiler:    " "${MPI_C_COMPILER}")
  message(STATUS "  MPI C flags:       " "${MPI_C_COMPILE_FLAGS}")
  message(STATUS "  MPI C include dir: " "${MPI_C_INCLUDE_PATH}")
  message(STATUS "  MPI C libraries:   " "${MPI_C_LIBRARIES}")
  message(STATUS "  MPI C flags:       " "${MPI_C_LINK_FLAGS}")
  message(STATUS "------------------------------------------------")
endif(WITH_MPI)

# if(APPLE) # message(STATUS "Trying to find Accelerate library") #
# find_library(BLAS_LIBRARIES Accelerate HINTS
# /System/Library/Frameworks/Accelerate.framework) #
# MARK_AS_ADVANCED(BLAS_LIBRARIES) # set(BLAS_LINKER_FLAGS
# "${CMAKE_LINKER_FLAGS} -framework accelerate") # find_path(LAPACK_INCLUDE_DIRS
# #   NAMES clapack.h #   PATHS ${BLAS_INCLUDE_DIRS} #   NO_DEFAULT_PATH # ) #
# set(LAPACK_HEADER ${LAPACK_INCLUDE_DIRS}/clapack.h) # set(LAPACK_SUFFIX "_") #
# set(LAPACK_PREFIX)

# # message(STATUS "Using Accelerate framework on macOS") message("using the
# local blas library") # option(BLAS_VERBOSE "Print some additional information
# during BLAS libraries detection" ON) set(BLAS_DIR "./mathlibs/src/blas/")
# set(BLA_VENDOR "Generic")

# set(ACCELERATE_LIBS "-framework Accelerate") BLAS_Accelerate_LIBRARY:FILEPATH=
# /System/Library/Frameworks/Accelerate.framework

# else()
if(WITH_MKL)
  # Find MKL
  find_package(MKL REQUIRED)
  set(HAVE_MKL TRUE)
  mark_as_advanced(HAVE_MKL)
  # Set BLAS and LAPACK to point to MKL
  set(BLAS_LIBRARIES
      ${MKL_BLAS_LIBRARIES}
      CACHE FILEPATH "")
  set(LAPACK_LIBRARIES
      ${MKL_LAPACK_LIBRARIES}
      CACHE FILEPATH "")
  if(MKL_SCALAPACK_LIBRARIES_FOUND)
    set(SCALAPACK_LIBRARIES ${MKL_SCALAPACK_LIBRARIES})
    # Versions >11.1 of MKL include Cluster PARDISO direct solver
    if(MKL_CPARDISO_FOUND)
      set(HAVE_CPARDISO TRUE)
      mark_as_advanced(HAVE_CPARDISO)
    endif()
  endif()
  include_directories(${MKL_INCLUDE_DIR})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_Fortran_FLAGS}")
else()
  if(NOT BLAS_LIBRARIES OR NOT LAPACK_LIBRARIES)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    # BLAS and LAPACK libraries may be needed for testing
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${BLAS_LIBRARIES}
                                 ${LAPACK_LIBRARIES})
    # If BLAS/LAPACK not explicitly set, attempt to find MKL BLAS/LAPACK
    # FIND_PACKAGE(MKL QUIET) IF (MKL_BLAS_LIBRARIES_FOUND AND
    # MKL_LAPACK_LIBRARIES_FOUND) # Set BLAS and LAPACK to point to MKL, in
    # order to work around # CMake FindBLAS finding wrong core library for GNU
    # Fortran SET(BLAS_LIBRARIES ${MKL_BLAS_LIBRARIES} CACHE FILEPATH "")
    # SET(LAPACK_LIBRARIES ${MKL_LAPACK_LIBRARIES} CACHE FILEPATH "")
    # INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR}) SET(CMAKE_Fortran_FLAGS
    # "${CMAKE_Fortran_FLAGS} ${MKL_Fortran_FLAGS}") ENDIF()
  endif()
endif()
# endif(APPLE)

message(STATUS "------------------------------------------------")
message(STATUS "  BLAS library:   " "${BLAS_LIBRARIES}")
message(STATUS "  LAPACK library: " "${LAPACK_LIBRARIES}")
message(STATUS "------------------------------------------------")

# MARK_AS_ADVANCED(WITH_Trilinos)
if(WITH_Trilinos)
  find_package(ML REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/ML ${TRILINOS_PATH})
  find_package(Trilinos REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos
               ${TRILINOS_PATH})
  find_package(Belos REQUIRED PATHS ${TRILINOS_PATH}/lib/cmake/Belos
               ${TRILINOS_PATH})
  if(ML_DIR)
    add_definitions(-DHAVE_TRILINOS)
  endif()
  message(STATUS "  Trilinos:           " "${Trilinos_FOUND}")
  message(STATUS "  Trilinos_DIR:       " "${Trilinos_DIR}")
  message(STATUS "  Trilinos_LIBRARIES: " "${Trilinos_LIBRARIES}")
  message(STATUS "------------------------------------------------")
endif()

if(WITH_Zoltan)
  # ZOLTAN Remeshing
  message(STATUS "------------------------------------------------")

  set(USE_SYSTEM_ZOLTAN
      FALSE
      CACHE BOOL "Use system Zoltan library")
  if(NOT ${USE_SYSTEM_ZOLTAN})
    message(STATUS "Building Zoltan mesh repartitioner")
    add_subdirectory(contrib/Zoltan_v3.83)
    set(HAVE_ZOLTAN
        TRUE
        CACHE BOOL "Has Zoltan tools for mesh (re)partitioning")
    mark_as_advanced(HAVE_ZOLTAN)
    message(STATUS "  Zoltan_LIB:      " "${ZOLTAN_LIBRARY}")
    message(STATUS "------------------------------------------------")
  else()
    find_package(Zoltan)

    if(ZOLTAN_FOUND)

      set(HAVE_ZOLTAN
          TRUE
          CACHE BOOL "Has Zoltan tools for mesh (re)partitioning")
      get_filename_component(ZOLTAN_LIBDIR ${ZOLTAN_LIBRARY} DIRECTORY)

      mark_as_advanced(HAVE_ZOLTAN)
      message(STATUS "  Zoltan:           " "${ZOLTAN_FOUND}")
      message(STATUS "  Zoltan_INC:       " "${ZOLTAN_INCLUDE_DIR}")
      message(STATUS "  Zoltan_LIB:      " "${ZOLTAN_LIBRARY}")
      message(STATUS "  Zoltan_LIBDIR:      " "${ZOLTAN_LIBDIR}")

      include_directories(${ZOLTAN_INCLUDE_DIR})
      message(STATUS "------------------------------------------------")

    else()
      message(STATUS "  Library not found: >ZOLTAN_FOUND< ")
      message(
        STATUS
          "    Missing: >ZOLTAN_INCLUDE_DIR< , >ZOLTAN_LIBRARY< for mesh (re)partitioning"
      )
    endif(ZOLTAN_FOUND)

    message(STATUS "------------------------------------------------")
  endif()
endif(WITH_Zoltan)

if(WITH_MMG)

  # MMG Remeshing
  message(STATUS "------------------------------------------------")
  message(STATUS "Mesh adaptation 2D/3D looking for [Mmg] tools ")

  find_package(MMG)
  # INCLUDE("cmake/Modules/FindMMG.cmake")
  # INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR}) TARGET_LINK_LIBRARIES(
  # ${YOUR_TARGET} ${MMG_LIBRARY})

  if(MMG_FOUND)

    set(HAVE_MMG
        TRUE
        CACHE BOOL "Has Mmg tools for Elmer/Ice Mesh adaptation")
    get_filename_component(MMG_LIBDIR ${MMG_LIBRARY} DIRECTORY)

    mark_as_advanced(HAVE_MMG)
    message(STATUS "  Mmg:           " "${MMG_FOUND}")
    message(STATUS "  Mmg_INC:       " "${MMG_INCLUDE_DIR}")
    message(STATUS "  Mmg_LIB:      " "${MMG_LIBRARY}")
    message(STATUS "  Mmg_LIBDIR:      " "${MMG_LIBDIR}")
    message(STATUS "Compile MMG2DSolver/MMG3DSolver")
    add_definitions(-DHAVE_MMG)

    include_directories(${MMG_INCLUDE_DIR})
    message(STATUS "------------------------------------------------")

  else()
    message(STATUS "  Library not found: >MMG_FOUND< ")
    message(
      STATUS
        "    Missing: >MMG_INCLUDE_DIR< , >MMG_LIBRARY<, to compile MMG3DSolver"
    )
  endif(MMG_FOUND)

  message(STATUS "------------------------------------------------")
endif()

if(WITH_Mumps)
  find_package(Mumps REQUIRED)
  if(Mumps_FOUND)
    set(HAVE_MUMPS TRUE)
    mark_as_advanced(HAVE_MUMPS)
    # LINK_DIRECTORIES(${MUMPS_LIBRARIES})
  endif()
  message(STATUS "  Mumps:             " "${Mumps_FOUND}")
  message(STATUS "  Mumps include:     " "${Mumps_INCLUDE_DIR}")
  message(STATUS "  Mumps libraries:   " "${Mumps_LIBRARIES}")
  message(STATUS "------------------------------------------------")
endif()

# Hypre
if(WITH_Hypre)
  find_package(
    Hypre REQUIRED
    OPTIONAL_COMPONENTS
      Euclid
      parcsr_ls
      ParaSails
      IJ_mv
      FEI
      core
      QUIET)
  if(Hypre_FOUND)
    set(HAVE_HYPRE TRUE)
    mark_as_advanced(HAVE_HYPRE)
    include_directories(${Hypre_INCLUDE_DIR})
    # LINK_DIRECTORIES(${Hypre_LIBRARIES})
  endif()
  message(STATUS "  Hypre:             " "${Hypre_FOUND}")
  message(STATUS "  Hypre include:     " "${Hypre_INCLUDE_DIR}")
  message(STATUS "  Hypre libraries:   " "${Hypre_LIBRARIES}")
  message(STATUS "------------------------------------------------")
endif()

# Permon
if(WITH_FETI4I)
  find_package(FETI4I REQUIRED)
  if(FETI4I_FOUND)
    set(HAVE_FETI4I TRUE)
    mark_as_advanced(HAVE_FETI4I)
    add_definitions(-DHAVE_FETI4I)
  endif()
endif()

# Check if Fortran compiler supports procedure pointer
include(testProcedurePointer)
if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_PROCEDUREPOINTER)
  message(
    FATAL_ERROR
      "Fortran compiler does not seem to support the PROCEDURE statement.")
endif()

# Check if Fortran compiler supports contiguous keyword
include(testContiguous)
if(USE_CONTIGUOUS AND CMAKE_Fortran_COMPILER_SUPPORTS_CONTIGUOUS)
  add_definitions(-DCONTIG=,CONTIGUOUS)
else()
  add_definitions(-DCONTIG=)
endif()

# Check if this is mingw toolchain.
if(MINGW)
  add_definitions(-DMINGW32)
endif()

if(WIN32)
  add_definitions(-DWIN32)
endif()

# Check if Fortran compiler supports execute_command_line
include(testExecutecommand)
if(CMAKE_Fortran_COMPILER_SUPPORTS_EXECUTECOMMANDLINE)
  add_definitions(-DHAVE_EXECUTECOMMANDLINE)
endif()

add_definitions(-DUSE_ISO_C_BINDINGS)
set(FC_STDCALLBULL " ")
set(ELMER_LINKTYP 1)
# SET(FC_FUNC FC_GLOBAL) SET(FC_FUNC_ FC_GLOBAL_)

# Add support for arpack
add_definitions(-DUSE_ARPACK)

# INCLUDE(FortranCInterface) FortranCInterface_HEADER(FCMangle.h MACRO_NAMESPACE
# "FC_") FortranCInterface_VERIFY(CXX)

# Set up rpaths to point ELMER_INSTALL_LIB_DIR
mark_as_advanced(ELMER_INSTALL_LIB_DIR)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(ELMERSOLVER_RPATH_STRING "\$ORIGIN/../${ELMER_INSTALL_LIB_DIR}")
  set(ELMERLIB_RPATH_STRING "\$ORIGIN/")
elseif(APPLE)
  set(ELMERSOLVER_RPATH_STRING "@loader_path/../${ELMER_INSTALL_LIB_DIR}")
  set(ELMERLIB_RPATH_STRING "@loader_path/")
endif()

# Uncomment these to use rpath with linked libraries
# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) SET(CMAKE_INSTALL_RPATH
# "${CMAKE_INSTALL_PREFIX}/lib") Remove system directories (from the CMake
# guide) LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
# "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir) IF("${isSystemDir}" STREQUAL "-1")
# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib") ENDIF("${isSystemDir}"
# STREQUAL "-1")

if(WITH_ElmerIce)
  message(STATUS "Adding optional package ElmerIce")
  add_subdirectory(elmerice)
endif(WITH_ElmerIce)

find_package(METIS REQUIRED)

add_subdirectory(matc)
add_subdirectory(umfpack)
add_subdirectory(fhutiter)
add_subdirectory(meshgen2d)
add_subdirectory(fem)
add_subdirectory(mathlibs)
add_subdirectory(elmergrid)
add_subdirectory(license_texts)

if(WITH_ELMERGUI)
  message(STATUS "  Building ElmerGUI")
  message(STATUS "------------------------------------------------")
  add_subdirectory(ElmerGUI)
endif(WITH_ELMERGUI)

if(WITH_ELMERGUITESTER)
  add_subdirectory(ElmerGUItester)
  message(STATUS "  Building ElmerGUI tester")
  message(STATUS "------------------------------------------------")
endif(WITH_ELMERGUITESTER)

if(WITH_ELMERGUILOGGER)
  add_subdirectory(ElmerGUIlogger)
  message(STATUS "  Building ElmerGUI logger")
  message(STATUS "------------------------------------------------")
endif(WITH_ELMERGUILOGGER)

# INCLUDE(FeatureSummary) FEATURE_SUMMARY(WHAT ENABLED_FEATURES
# INCLUDE_QUIET_PACKAGES DESCRIPTION "Enabled Features:" VAR
# enabledFeaturesText) MESSAGE(STATUS "${enabledFeaturesText}")

#
if(CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPCASE)
  mark_as_advanced(CMAKE_BUILD_TYPE_UPCASE)
endif(CMAKE_BUILD_TYPE)

message(STATUS "  Fortran compiler:        " "${CMAKE_Fortran_COMPILER}")
message(
  STATUS
    "  Fortran flags:           "
    "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
message(STATUS "------------------------------------------------")
message(STATUS "  C compiler:              " "${CMAKE_C_COMPILER}")
message(STATUS "  C flags:                 "
               "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
message(STATUS "------------------------------------------------")
message(STATUS "  CXX compiler:            " "${CMAKE_CXX_COMPILER}")
message(
  STATUS "  CXX flags:               "
         "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
message(STATUS "------------------------------------------------")

if(WITH_CONTRIB)
  message(STATUS "  Building contributed solvers")
  message(STATUS "------------------------------------------------")
endif(WITH_CONTRIB)

if(WITH_ELMERPOST)
  mark_as_advanced(WITH_ELMERPOST)
  message(WARNING "  Building ElmerPost")
  message(STATUS "------------------------------------------------")
  add_subdirectory(post)
endif(WITH_ELMERPOST)

if(HAVE_LUA)
  message(STATUS "  Building Lua")
  message(STATUS "------------------------------------------------")
endif(HAVE_LUA)

# Packaging
if(NOT BYPASS_CPACK)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cpack/ElmerCPack.cmake)
endif(NOT BYPASS_CPACK)
