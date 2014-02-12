## Makefile.am for hutiter/src

INCLUDES = \
        -I$(top_builddir)/include       \
        -I$(top_srcdir)/include         

COMMON_SRCS =  \
	eio-config.h \
	EIOConstraintAgent.cpp \
	EIODualMeshAgent.cpp \
	EIOGeometryAgent.cpp \
	EIOMeshAgent.cpp \
	EIOModelDataAgent.cpp \
	EIOModelManager.cpp \
	EIOPartWriter.cpp \
	EIOSolverAgent.cpp 

libeioc_a_SOURCES = $(COMMON_SRCS) eio_api_c.cpp
libeiof_a_SOURCES = $(COMMON_SRCS) eio_api_f.cpp

# if USE_SHARED_LIBS
# noinst_LIBRARIES = libeioc.a libeiof.a
# all-local: libeioc.$(SHLEXT) libeiof.$(SHLEXT)

# libeioc.$(SHLEXT): libeioc.a
#	$(SH_LD2) $(SH_LDFLAGS) $(B64FLAGS) $(LDFLAGS) -o $@ libeioc.a

# libeiof.$(SHLEXT): libeiof.a
# 	$(SH_LD2) $(SH_LDFLAGS) $(B64FLAGS) $(LDFLAGS) -o $@ libeiof.a
# else
lib_LIBRARIES = libeioc.a libeiof.a
# endif
