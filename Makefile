# Makefile to build the whole Elmer project with MPI, Scalapack, MUMPS and SuiteSparse

# Installation directory
IDIR = /home/mbycklin/code/elmerfem_git/bin/

# PARALLEL COMPILATION VARIABLES

# Parallel compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif90

# Sequential compilers
# export CC=gcc
# export CXX=gcc
# export FC=gfortran
# export F77=gfortran

# Parallel compilation variables
# export MPI_HOME=/opt/openmpi-intel
# export MPI_INCLUDE=-I/opt/openmpi-intel/include
# export MPI_HOME=/usr/lib64/openmpi
# export MPI_INCLUDE=-I/usr/include/openmpi-x86_64/
export MPI_HOME=/opt/mpich/3.0.4/
export MPI_INCLUDE=$(MPI_HOME)/include

# General compilation flags
OPT_FLAGS=-O3 -g -m64 -fopenmp -ftree-vectorize -funroll-loops
# OPT_FLAGS=-O0 -g -m64 -fopenmp -fbacktrace -Warray-bounds -Wuninitialized -Wall
# OPT_FLAGS=-O0 -g -m64 -fopenmp -fbacktrace
# OPT_FLAGS=-O3 -m64 -openmp

# BLAS configure flags
BLASCFLAGS=--with-blas=$(HOME)/code/GotoBLAS2/libgoto2_penryn-r1.13.so
LAPACKCFLAGS=--with-lapack=$(HOME)/code/GotoBLAS2/libgoto2_penryn-r1.13.so

# Parallel Elmer configure flags
ELMERCFLAGS=--with-mpi=yes --with-mpi-dir=$(MPI_HOME) --with-64bits=yes
ELMERTARGETS=ElmerGrid Mesh2D matc elmerf90 elmerld ElmerSolver \
             GebhardtFactors ViewFactors elmerf90-nosh ElmerMesh2D \
             ElmerSolver_mpi SC2Elmer

export ELMER_HOME=$(IDIR)

# MKL related flags
# MKLLIB=mkl
ifdef MKLLIB
MKLFLAGS=-DHAVE_MKL -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
MKLFCLAGS=
MKLLDFLAGS=-L$(MKLROOT)/lib/intel64 -lmkl_blas95_lp64 \
           -lmkl_lapack95_lp64 -lmkl_scalapack_lp64 \
           -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 \
           -openmp -lpthread -lm
ELMERCFLAGS=--with-mpi=yes --with-mpi-dir=$(MPI_HOME) --with-64bits=yes \
            --with-blas=-mkl \
            --with-lapack=-mkl
endif

# Scalapack related flags (uncomment SCALAPACKLIB to compile Scalapack)
ifndef MKLLIB
# SCALAPACKLIB=scalapack
# SCALAPACKPATH=/home/mbycklin/code/elmerfem_git/scalapack
# SCALAPACKLDFLAGS=
endif

# Mumps related flags (uncomment MUMPSLIB to compile Mumps)
# MUMPSLIB=mumps
# MUMPSPATH=/home/mbycklin/code/elmerfem_git/mumps
# MUMPSCFLAGS=-DHAVE_MUMPS
# MUMPSFCFLAGS=-I$(MUMPSPATH)/include
ifdef MKLLIB
MUMPSLDFLAGS=-L$(MUMPSPATH)/lib -ldmumps -lmumps_common \
             -lpord -mkl \
             -L$(MPI_HOME)/lib -lmpi_f77 -lmpi
endif
ifdef SCALAPACKLIB
MUMPSLDFLAGS=-L$(MUMPSPATH)/lib -ldmumps -lmumps_common \
             -lpord -L$(SCALAPACKPATH)/ -lscalapack \
             -L$(MPI_HOME)/lib -lmpi_f77 -lmpi -lblas
endif

# SuiteSparse related flags (uncomment SUITESPARSELIB to compile SuiteSparse)
# SUITESPARSELIB=suitesparse
ifdef SUITESPARSELIB
SUITESPARSEPATH=/home/mbycklin/code/elmerfem_git/suitesparse
SUITESPARSECFLAGS=-DHAVE_UMFPACK -DHAVE_SPQR -DHAVE_CHOLMOD
ifdef MKLLIB
SUITESPARSELDFLAGS=-L$(SUITESPARSEPATH)/UMFPACK/Lib -lumfpack \
                   -L$(SUITESPARSEPATH)/SPQR/Lib -lspqr \
                   -L$(SUITESPARSEPATH)/CHOLMOD/Lib -lcholmod \
                   -L$(SUITESPARSEPATH)/COLAMD/Lib -lcolamd \
                   -L$(SUITESPARSEPATH)/SuiteSparse_config -lsuitesparseconfig
else
SUITESPARSELDFLAGS=-L$(SUITESPARSEPATH)/UMFPACK/Lib -lumfpack \
                   -L$(SUITESPARSEPATH)/SPQR/Lib -lspqr \
                   -L$(SUITESPARSEPATH)/CHOLMOD/Lib -lcholmod \
                   -L$(SUITESPARSEPATH)/COLAMD/Lib -lcolamd \
                   -L$(SUITESPARSEPATH)/SuiteSparse_config -lsuitesparseconfig \
                   -llapack
endif
endif

# MatrixMarket IO library
# MMIOLIB=mmio
MMIOPATH=/home/mbycklin/code/elmerfem_git/mmio
ifdef MMIOLIB
MMIOCFLAGS=-DHAVE_MMIO
MMIOLDFLAGS=-L$(MMIOPATH)/lib -lmmio
endif

# Export final compilation flags
export CFLAGS=$(OPT_FLAGS) $(MUMPSCFLAGS) $(SUITESPARSECFLAGS) $(MKLFLAGS) $(MMIOCFLAGS)
export CPPFLAGS=$(OPT_FLAGS) $(MUMPSCFLAGS) $(SUITESPARSECFLAGS) $(MKLFLAGS) $(MMIOCFLAGS)
export FCPPFLAGS=$(OPT_FLAGS) $(MUMPSCFLAGS) $(SUITESPARSECFLAGS) $(MKLFLAGS) $(MMIOCFLAGS)
export CXXFLAGS=$(OPT_FLAGS) $(MUMPSCFLAGS) $(SUITESPARSECFLAGS) $(MKLFLAGS) $(MMIOCFLAGS)
export FCFLAGS=$(OPT_FLAGS) $(MUMPSFCFLAGS) $(SUITESPARSECFLAGS) $(MKLFLAGS) $(MKLFCLAGS) $(MMIOCFLAGS)
export F90FLAGS=$(OPT_FLAGS) $(MKLFLAGS) $(MKLFCLAGS) $(MMIOCFLAGS)
export F77FLAGS=$(OPT_FLAGS) $(MKLFLAGS) $(MKLFCLAGS) $(MMIOCFLAGS)
export FFLAGS=$(OPT_FLAGS) $(MKLFLAGS) $(MKLFCLAGS) $(MMIOCFLAGS)

# Final Elmer library flags
ELMERLDFLAGS=$(MKLLDFLAGS) $(MUMPSLDFLAGS) $(SCALAPACKLDFLAGS) $(SUITESPARSELDFLAGS) $(MMIOLDFLAGS)

# Library names (required)
# RLIBDIRS=matc umfpack mathlibs elmergrid meshgen2d eio hutiter umfpack
RLIBDIRS=matc umfpack mathlibs elmergrid meshgen2d eio hutiter
# Library names (optional)
OLIBDIRS=$(SCALAPACKLIB) $(MUMPSLIB) $(SUITESPARSELIB)

LIBDIRS=$(RLIBDIRS) $(OLIBDIRS)

# MAKE RULES
.PHONY: test fem/test %/build %/install inst-clean
.PRECIOUS: %/Makefile $(MUMPSLIB)/Makefile.inc $(SCALAPACKLIB)/SLmake.inc \
           $(SUITESPARSELIB)/SuiteSparse_config/SuiteSparse_config.mk

# ElmerSolver build rules
all: lib-install fem/install
build: lib-build fem/build
install: lib-install fem/install
clean: lib-clean fem/clean inst-clean
distclean: lib-distclean fem/distclean inst-clean 
test: fem/test

# Library build rules 
lib-build: $(addsuffix /build, $(LIBDIRS))
lib-install: $(addsuffix /install, $(LIBDIRS))
lib-clean: $(addsuffix /clean, $(LIBDIRS))
lib-distclean: $(addsuffix /distclean, $(LIBDIRS))

# Generic configure rule
%/Makefile:
	cd $(dir $@); ./configure --prefix=$(IDIR) $(ELMERCFLAGS)

# Generic build rule
%/build: %/Makefile
	cd $(dir $@); make

# Generic install rule
%/install: %/build
	cd $(dir $@); make install

# Generic clean rule
%/clean: %/Makefile
	-cd $(dir $@); make clean

# Generic distclean rule
%/distclean: %/Makefile
	cd $(dir $@); make clean; cd ..; rm $(dir $@)Makefile

# More specific rules

# SCALAPACK
ifdef SCALAPACKLIB
$(SCALAPACKLIB)/SLmake.inc:
	cp slmake.inc.default $(SCALAPACKLIB)/SLmake.inc

$(SCALAPACKLIB)/build: $(SCALAPACKLIB)/SLmake.inc
	cd $(SCALAPACKLIB); make lib;

$(SCALAPACKLIB)/install: $(SCALAPACKLIB)/SLmake.inc
	-cd $(SCALAPACKLIB); make lib

$(SCALAPACKLIB)/distclean: $(SCALAPACKLIB)/SLmake.inc
	-cd $(SCALAPACKLIB); make clean; rm SLmake.inc
endif

# MUMPS
ifdef MUMPSLIB
$(MUMPSLIB)/Makefile.inc:
	cp mmake.inc.default $(MUMPSLIB)/Makefile.inc

$(MUMPSLIB)/build: $(MUMPSLIB)/Makefile.inc
	cd $(MUMPSLIB); make d

$(MUMPSLIB)/install: $(MUMPSLIB)/Makefile.inc
	cd $(MUMPSLIB); make d

$(MUMPSLIB)/distclean: $(MUMPSLIB)/Makefile.inc
	-cd $(MUMPSLIB); make clean; rm Makefile.inc
endif

# SuiteSparse
ifdef SUITESPARSELIB
$(SUITESPARSELIB)/SuiteSparse_config/SuiteSparse_config.mk:
	cp ssmake.mk.default $(SUITESPARSELIB)/SuiteSparse_config/SuiteSparse_config.mk

$(SUITESPARSELIB)/build: $(SUITESPARSELIB)/SuiteSparse_config/SuiteSparse_config.mk
	cd $(SUITESPARSELIB); make library

$(SUITESPARSELIB)/install: $(SUITESPARSELIB)/SuiteSparse_config/SuiteSparse_config.mk $(SUITESPARSELIB)/build 
	cd $(SUITESPARSELIB); make install

$(SUITESPARSELIB)/clean: $(SUITESPARSELIB)/SuiteSparse_config/SuiteSparse_config.mk
	-cd $(SUITESPARSELIB); make clean

$(SUITESPARSELIB)/distclean: $(SUITESPARSELIB)/SuiteSparse_config/SuiteSparse_config.mk
	-cd $(SUITESPARSELIB); make purge; rm SuiteSparse_config/SuiteSparse_config.mk

endif

# MatrixMarket IO
ifdef MMIOLIB
$(MMIOLIB)/build:
	cd $(MMIOLIB); make lib

$(MMIOLIB)/install:
	cd $(MMIOLIB); make lib

$(MMIOLIB)/clean: 
	-cd $(MMIOLIB); make clean
endif

# fem (elmersolver) library
fem/Makefile:
	cd fem; ./configure --prefix=$(IDIR) $(ELMERCFLAGS) $(BLASCFLAGS) $(LAPACKCFLAGS) LDFLAGS="$(ELMERLDFLAGS)"

fem/build: fem/Makefile
	cd fem; make

fem/install: fem/build
	cd fem; make install

fem/test:
	cd fem/tests; make check

fem/clean:
	-cd fem; make clean

fem/distclean:
	-cd fem; make distclean

# Installation clean
inst-clean:
	rm -f $(addprefix $(IDIR)bin/,$(ELMERTARGETS)); \
	rm -f $(IDIR)include/*.h;
	rm -f $(IDIR)include/elmer/*.h; \
	rm -f $(IDIR)lib/*.a $(IDIR)lib/*.so; \
	rm -rf $(IDIR)share/elmersolver; \
