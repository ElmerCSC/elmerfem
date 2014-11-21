
########### next target ###############
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(visuals_STAT_SRCS
   arrows.c
   colscale.c
   contour_lines.c
   isosurface.c
   mesh.c
   particles.c
   sphere.c
   visual.c
   visual.h
)

add_library(visuals STATIC ${visuals_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

#noinst_LIBRARIES=libvisuals.a
#
#libvisuals_a_SOURCES=\
#	arrows.c \
#	colscale.c \
#	contour_lines.c \
#	isosurface.c \
#	mesh.c \
#	particles.c \
#	sphere.c \
#	visual.c \
#	visual.h 
#
