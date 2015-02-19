# Test for Grounding Line
#

PLATFORM := $(shell uname -s)

run:
	if [ '$(PLATFORM)' = 'Linux' ] ; \
	then \
           cp $(ELMER_HOME)/share/elmersolver/lib/FreeSurfaceSolver.so MyFreeSurfaceSolver.so ; \
	else \
	   cp $(ELMER_HOME)/share/elmersolver/lib/FreeSurfaceSolver.dylib MyFreeSurfaceSolver.dylib ; \
        fi ;
	$(ELMER_SOLVER) mismip.sif

clean:
	/bin/rm *.dat *.log *.dat.names GL_MISMIP.txt
	/bin/rm MyFreeSurfaceSolver.*
