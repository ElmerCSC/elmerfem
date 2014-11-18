include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/include )


########### next target ###############

SET(camera_STAT_SRCS
   camera.c
   glp.cc
   file.cc
   context.cc
)

add_library(camera STATIC ${camera_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

## 
## jv: merged libGLP with camera.
## 
#noinst_LIBRARIES=libcamera.a
#
#libcamera_a_SOURCES=camera.c glp.cc file.cc context.cc 
#noinst_HEADERS=camera.h glpfile.h glp.h
#
#INCLUDES = \
#        -I$(top_srcdir)/src \
#        -I$(top_srcdir)/src/include
#
