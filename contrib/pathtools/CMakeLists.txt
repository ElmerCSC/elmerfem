# The MinGW pathtools have originally been written by Ray Donnelly in 2014
# under a permissive license.
# They have been used since in MSYS2 to help with "relocation tasks".
# They are mainly written for and tested on Windows. They might not work
# correctly on other platforms.

# Create an object library so the pathtools can be used in different
# executables or libraries.
ADD_LIBRARY(pathtools OBJECT "pathtools.c")

TARGET_INCLUDE_DIRECTORIES(pathtools
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
