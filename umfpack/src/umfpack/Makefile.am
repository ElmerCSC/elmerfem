SUBDIRS=include

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------
include make.dep

EXTRA_DIST = $(UMFPACK_ORIG_SOURCES)
UMFPACK_SOURCES = umfpack_timer.c umfpack_tictoc.c

INCLUDES = \
	-I./include       \
	-I$(top_srcdir)/src/amd/include     \
	-I$(top_srcdir)/src/amd \
	-I.         

lib_LIBRARIES = libumfpack.a

libumfpack_a_SOURCES = $(UMFPACK_SOURCES)
nodist_libumfpack_a_SOURCES = $(UMFPACK_CPP_SOURCES)

clean-local:
	-rm -f $(UMFPACK_CPP_SOURCES)
