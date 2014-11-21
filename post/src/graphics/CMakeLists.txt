
########### next target ###############
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BUILD_DIR}/../)
GET_DIRECTORY_PROPERTY(local_inc_dirs INCLUDE_DIRECTORIES ) 
message(STATUS "graphics: inc_dirs = ${local_inc_dirs}")
SET(graphics_STAT_SRCS
   graphics.c
   graphics.h
   gutil1.c
   gutil2.c
)

add_library(graphics STATIC ${graphics_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

#noinst_LIBRARIES=libgraphics.a
#
#libgraphics_a_SOURCES=\
#	graphics.c \
#	graphics.h \
#	gutil1.c \
#	gutil2.c 
