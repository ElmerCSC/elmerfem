# Test for Grounding Line
#

PLATFORM := $(shell uname -s)

run:
	if [ '$(PLATFORM)' = 'Linux' ] ; \
	then \
           cp $(ELMER_HOME)/share/elmersolver/lib/FreeSurfaceSolver.so MyFreeSurfaceSolver ; \
	else \
	   cp $(ELMER_HOME)/share/elmersolver/lib/FreeSurfaceSolver.dylib MyFreeSurfaceSolver ; \
        fi ;
	$(ELMER_SOLVER) mismip.sif

clean:
	/bin/rm *.dat *.log *.dat.names
	/bin/rm MyFreeSurfaceSolver.*
