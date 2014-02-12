#-------------------------------------------------------------------------------
# AMD Makefile for compiling on Unix systems (for GNU make only)
#-------------------------------------------------------------------------------
SUBDIRS=include

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------
include make.dep

AMD_SOURCES = amd_internal.h 


INCLUDES = \
	-I./include       \
	-I.         

lib_LIBRARIES = libamd.a libamdf77.a

#-------------------------------------------------------------------------------
# Create the libamd.a library (C versions only)
#-------------------------------------------------------------------------------

libamd_a_SOURCES = $(AMD_SOURCES)
nodist_libamd_a_SOURCES = $(AMD_CPP_SOURCES)

#-------------------------------------------------------------------------------
# compile the Fortran versions and the libamdf77.a library
#-------------------------------------------------------------------------------

AMDF77_SOURCES = amd.f amdbar.f
libamdf77_a_SOURCES = $(AMDF77_SOURCES)

clean-local:
	-rm -f $(AMD_CPP_SOURCES)

EXTRA_DIST= \
	amd_aat.c \
	amd_1.c \
	amd_2.c \
	amd_dump.c \
	amd_postorder.c \
	amd_post_tree.c \
	amd_defaults.c \
	amd_order.c \
	amd_control.c \
	amd_info.c \
	amd_valid.c \
	amd_preprocess.c
