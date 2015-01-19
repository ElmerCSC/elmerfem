
SHELL = /bin/sh
RM=/bin/rm
CP=/bin/cp

TARGETDIR = EXECUTABLES/.
ELMERDIR = $(ELMER_HOME)/share/elmersolver/lib/.
ELMERF90 = elmerf90
libs = \
$(TARGETDIR)/FreeSurfaceSolver1\
$(TARGETDIR)/FreeSurfaceSolver2\

run: $(libs)
	$(ELMER_GRID) 1 2 mesh.grd
	$(ELMER_SOLVER) damage.sif


clean: 
	$(RM) -f $(libs)
	$(RM) *.dat *.log *.dat.names
	$(RM) Damage.*
	$(RM) -r mesh

$(TARGETDIR)/FreeSurfaceSolver1: $(ELMERDIR)/FreeSurfaceSolver.* 
	$(CP) $(ELMERDIR)/FreeSurfaceSolver.so $(TARGETDIR)/FreeSurfaceSolver1

$(TARGETDIR)/FreeSurfaceSolver2: $(ELMERDIR)/FreeSurfaceSolver.* 
	$(CP) $(ELMERDIR)/FreeSurfaceSolver.so $(TARGETDIR)/FreeSurfaceSolver2
