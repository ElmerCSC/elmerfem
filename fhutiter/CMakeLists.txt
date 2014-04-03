# For the installer process we need a single directory
# where we put all module files
SET(CMAKE_Fortran_MODULE_DIRECTORY
  ${PROJECT_BINARY_DIR}/fmodules CACHE PATH "Directory for Fortran modules")

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/fhutiter/src")

ADD_SUBDIRECTORY(src)
