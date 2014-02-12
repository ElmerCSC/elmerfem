#### Makefile to compile the solver Scattered2DDataInterpolator
##### This solver require external libraries:
	#####   netcdf
	#####	csa: http://code.google.com/p/csa-c/
	#####	nn: http://code.google.com/p/nn-c/
####### These libraries must be compiled by the user and this Makefile edited/updated


####### EDIT the floowing flags according to your installation ###########
###### C compilert and FLAG
CC=mpicc
CCFLAG=-fPIC

# PATH to the HOME directory of NN and CSA header files and libraries
NN_HOME=/Users/fgillet/Programmes/INTERPOLATION/nn
CS_HOME=/Users/fgillet/Programmes/INTERPOLATION/csa
NNLIB=$(NN_HOME)/libnn.a
CSLIB=$(CS_HOME)/libcsa.a

# NETCDF INCLUDE and LIBS Flags
NETCDFHOME="/Users/fgillet/Programmes/install"
NETCDFINC="-I$(NETCDFHOME)/include"
NETCDFLIBS="-L$(NETCDFHOME)/lib -lnetcdff -lnetcdf"




##########################################################################
ELMERICE_HOME ?= $(ELMER_HOME)/share/elmersolver
##############

all: Scattered2DDataInterpolator

install: Scattered2DDataInterpolator 
	cp -f Scattered2DDataInterpolator $(ELMERICE_HOME)/lib/.

Scattered2DDataInterpolator: Scattered2DDataInterpolator.f90 csa_interpolate_points.o Scattered2D_FInterface.o
	elmerf90 $(NETCDFINC) $(NNLIB) $(CSLIB)  $^ -o $@ $(NETCDFLIBS)

csa_interpolate_points.o : csa_interpolate_points.c
	$(CC) $(CCFLAG) -I$(CS_HOME) -c $^ -o $@

Scattered2D_FInterface.o : Scattered2D_FInterface.f90
	elmerf90-nosh -c $< -o $@

clean:
	rm -f Scattered2DDataInterpolator
	rm -f *.o
	rm -f *.mod
