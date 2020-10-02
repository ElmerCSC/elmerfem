SrcDir=./src
ExeDir=./Executables


_SOURCES=USF_Init.f90  GradientValidation.f90
SOURCES = $(patsubst %,$(SrcDir)/%,$(_SOURCES))

PLATFORM := $(shell uname -s)
ifeq  ('$(PLATFORM)', 'Darwin') 
   _SOLVERS=$(_SOURCES:.f90=.dylib)
else
   _SOLVERS=$(_SOURCES:.f90=.so)
endif
SOLVERS=$(patsubst %,$(ExeDir)/%,$(_SOLVERS)) 

all:$(SOLVERS)

Mesh: 
	ElmerGrid  1 2 mesh2D -metis $(np) 4

$(ExeDir)/%.dylib: $(SrcDir)/%.f90
	elmerf90 $< -o $@

$(ExeDir)/%.so: $(SrcDir)/%.f90
	elmerf90 $< -o $@


clean:
	rm -f *~
	rm -f *.dat* 
	rm -f *.out
	rm -f Executables/*
	rm -f -r mesh2D
	rm -f ELMERSOLVER_STARTINFO
