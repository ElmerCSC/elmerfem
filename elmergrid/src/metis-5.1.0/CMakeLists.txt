cmake_minimum_required(VERSION 2.8)
project(METIS)

#set(GKLIB_PATH "GKlib" CACHE PATH "path to GKlib")
set(GKLIB_PATH ${PROJECT_SOURCE_DIR}/GKlib CACHE PATH "GKlib source path")
set(METIS_SHARED FALSE CACHE BOOL "build a shared library")

# if(MSVC)
#   set(METIS_INSTALL FALSE)
# else()
#   set(METIS_INSTALL TRUE)
# endif()

# Configure libmetis library.
if(METIS_SHARED)
  set(METIS_LIBRARY_TYPE SHARED)
else()
  set(METIS_LIBRARY_TYPE STATIC)
endif(METIS_SHARED)

include(${GKLIB_PATH}/GKlibSystem.cmake)
# Add include directories.
include_directories(BEFORE ${GKLIB_PATH})
include_directories(BEFORE include)
# Recursively look for CMakeLists.txt in subdirs.
add_subdirectory("include")
add_subdirectory("libmetis")
#add_subdirectory("programs")
