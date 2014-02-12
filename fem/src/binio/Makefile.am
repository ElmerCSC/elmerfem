SUBDIRS = test

#AM_CFLAGS=-Wall -ansi -pedantic -O3
AM_FCFLAGS=@INCLUDE_MODULE_FLAG@.  @INCLUDE_MODULE_FLAG@..

noinst_LIBRARIES = libbinio.a
libbinio_a_SOURCES = binio.c kinds.f90 biniomod.f90

biniomod.f90: kinds.o

ELMER_SOLVER_DATADIR=$(prefix)/share/elmersolver
install-exec-local:
	mkdir -p $(DESTDIR)$(ELMER_SOLVER_DATADIR)/include
	$(INSTALL) *.mod $(DESTDIR)$(ELMER_SOLVER_DATADIR)/include

CLEANFILES = *.mod
