OBJ_EXT=.$(OBJEXT)
SHL_EXT=.$(SHLEXT)
LIB_EXT=.$(LIBEXT)
EXE_EXT=$(EXEEXT)

include make.bin

ELMER_TARGETS=$(BINARIES)

all-local: $(ELMER_TARGETS)

SUFFIXES: .f90 .so .src 

AM_FCFLAGS=@INCLUDE_MODULE_FLAG@../ @INCLUDE_MODULE_FLAG@./

if USE_CPP
# Use basic CPP
.src.f90:
	$(CPP) $(FCPPFLAGS) $< > $*.f90

.f90.$(OBJEXT):
	$(FC) $(AM_FCFLAGS)  $(FCFLAGS) -c $<
else
# Use fortran built-in cpp
.src.f90:
	$(CP) $< $*.f90

.f90.$(OBJEXT):
	$(FC) $(AM_FCFLAGS)  $(FCFLAGS) $(FORTRAN_CPP_FLAG) $(FCPPFLAGS) -c $<
endif

if USE_WINDOWS_COMPILER

.$(OBJEXT).$(SHLEXT):
	$(SH_LD) $(SH_LDFLAGS) $(SH_LINKING_TO_FLAGS) $(B64FLAGS) -o $@ $*$(OBJ_EXT) $(EXTRA_LIBS) -L.. -lelmersolver

EliminatePeriodic.$(SHLEXT): EliminateDirichlet$(OBJ_EXT) EliminatePeriodic$(OBJ_EXT)
	$(SH_LD) $(SH_LDFLAGS) $(SH_LINKING_TO_FLAGS) $(B64FLAGS) -o $@ $*$(OBJ_EXT) $(EXTRA_LIBS) EliminateDirichlet$(OBJ_EXT) -L.. -lelmersolver


else

.$(OBJEXT).$(SHLEXT): 
	$(SH_LD) $(SH_LDFLAGS) $(SH_LINKING_TO_FLAGS) $(B64FLAGS) -o $@ $*$(OBJ_EXT) $(EXTRA_LIBS)

endif


clean-local:
	$(RM) $(ALL) *.o *.mod *.f90 *.F *.fpp *~ *.a *.so *.dll \#*


ELMER_SOLVER_DATADIR=$(prefix)/share/elmersolver
install-data-local:
	@echo "Installing solvers to ${ELMER_SOLVER_DATADIR}/lib"

	mkdir -p $(DESTDIR)$(ELMER_SOLVER_DATADIR)/lib
	$(CP) *$(SHL_EXT) $(DESTDIR)$(ELMER_SOLVER_DATADIR)/lib

	@echo "Installing *.mod to ${ELMER_SOLVER_DATADIR}/include"
	mkdir -p $(DESTDIR)$(ELMER_SOLVER_DATADIR)/include
	$(CP) *.mod $(DESTDIR)$(ELMER_SOLVER_DATADIR)/include
