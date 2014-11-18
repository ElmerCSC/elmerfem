
########### next target ###############

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
