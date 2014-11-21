
########### next target ###############
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(tk_STAT_SRCS
   cursor.c
   event.c
   font.c
   getset.c
   image.c
   private.h
   window.c
   shapes.c
)

add_library(tk STATIC ${tk_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

#
#if NOTMINGW32
#
#noinst_LIBRARIES=libtk.a
#
#libtk_a_SOURCES=\
#	cursor.c \
#	event.c \
#	font.c \
#	getset.c \
#	image.c \
#	private.h \
#	window.c \
#	shapes.c 
#
#noinst_HEADERS=tk.h private.h
#
#endif
