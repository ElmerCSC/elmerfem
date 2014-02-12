SUBDIRS=view3d viewaxis binio . modules

OBJ_EXT=.$(OBJEXT)
SHL_EXT=.$(SHLEXT)
LIB_EXT=.$(LIBEXT)
EXE_EXT=$(EXEEXT)

include make.lib
include make.dep

NORMAL_TARGETS=mpif libelmersolver$(SHL_EXT) ElmerSolver$(EXE_EXT) ViewFactors$(EXE_EXT) GebhardtFactors$(EXE_EXT)
MPI_TARGETS=ElmerSolver_mpi$(EXE_EXT)


if USE_MPI
ELMER_TARGETS=$(NORMAL_TARGETS) $(MPI_TARGETS)
else
ELMER_TARGETS=$(NORMAL_TARGETS)
endif

all-local: $(ELMER_TARGETS)

SUFFIXES: .f90 .so .src .cpp

if USE_CPP
# Use basic CPP
.src.f90:
	$(CPP) $(FCPPFLAGS) $< > $*.f90

.f90.$(OBJEXT):
	$(FC) $(FCFLAGS) -c $<
else
# Use fortran built-in cpp
.src.f90:
	$(CP) $< $*.f90

.f90.$(OBJEXT):
	$(FC) $(FCFLAGS) $(FORTRAN_CPP_FLAG) $(FCPPFLAGS) -c $<
endif

.cxx.$(OBJEXT):
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

nofdeps:
	sed -e "s/^\([^:^#]*\...OBJEXT.: \)/#\1/g" Makefile > Makefile.nodeps
	mv Makefile.nodeps Makefile

# remove cross dependencies from fortran files
fdeps:
	sed -e 's/^#\([^:]*\...OBJEXT.: \)/\1/g' Makefile > Makefile.deps
	mv Makefile.deps Makefile


if USE_WINDOWS_COMPILER

.$(OBJEXT).$(SHLEXT): 
	$(SH_LD) $(SH_LDFLAGS) $(SH_LINKING_TO_FLAGS) $(B64FLAGS) -o $@ $*$(OBJ_EXT) $(EXTRA_LIBS) -L. -lelmersolver

EliminatePeriodic.$(SHLEXT): EliminateDirichlet$(OBJ_EXT) EliminatePeriodic$(OBJ_EXT)
	$(SH_LD) $(SH_LDFLAGS) $(SH_LINKING_TO_FLAGS) $(B64FLAGS) -o $@ $*$(OBJ_EXT) $(EXTRA_LIBS) EliminateDirichlet$(OBJ_EXT) -L. -lelmersolver

if USE_MPI
libelmersolver$(SHL_EXT): $(SOLVEROBJS) binio/libbinio.a
	$(RM) $@
	$(SH_LD2) $(SH_LDFLAGS) $(B64FLAGS) $(LDFLAGS) -o $@ $(SOLVEROBJS) $(SOLVER_LIBS) $(MPI_LIBS) -L. -Lbinio -lbinio
else
libelmersolver$(SHL_EXT): $(SOLVEROBJS) mpi_stubs$(OBJ_EXT) binio/libbinio.a
	$(RM) $@
	$(SH_LD2) $(SH_LDFLAGS) $(B64FLAGS) $(LDFLAGS) -o $@ $(SOLVEROBJS) $(SOLVER_LIBS) mpi_stubs$(OBJ_EXT) -L. -Lbinio -lbinio
endif

else

libelmersolver$(SHL_EXT): $(SOLVEROBJS) binio/libbinio.a
	$(RM) $@
if MAC_OSX
	$(SH_LD2) $(RPATH_ELMER) $(SH_LDFLAGS) $(B64FLAGS) $(LDFLAGS) -Wl,-dylib_install_name,libelmersolver-$(VERSION)$(SHL_EXT) -o libelmersolver-$(VERSION)$(SHL_EXT) $(SOLVEROBJS) $(SOLVER_LIBS) -L. -Lbinio -lbinio
else
if IBM_AIX
	$(SH_LD2) $(RPATH_ELMER) $(SH_LDFLAGS) $(B64FLAGS) $(LDFLAGS) -o libelmersolver-$(VERSION)$(SHL_EXT) $(SOLVEROBJS) $(SOLVER_LIBS) -L. -Lbinio -lbinio
else
	$(SH_LD2) $(RPATH_ELMER) $(SH_LDFLAGS) $(B64FLAGS) $(LDFLAGS)  -Wl,-soname,libelmersolver-$(VERSION)$(SHL_EXT) -o libelmersolver-$(VERSION)$(SHL_EXT) $(SOLVEROBJS) $(SOLVER_LIBS) -L. -Lbinio -lbinio
endif
endif
	$(LN) -s libelmersolver-$(VERSION)$(SHL_EXT) $@


.$(OBJEXT).$(SHLEXT): 
	$(SH_LD) $(SH_LDFLAGS) $(SH_LINKING_TO_FLAGS) $(B64FLAGS) -o $@ $*$(OBJ_EXT) $(EXTRA_LIBS)

endif


$(SOLVEROBJS): mpif

if USE_MPI
# the only way I got ibm MPI working was with mpxlf90...
if IBM_AIX
ElmerSolver_mpi$(EXE_EXT): Solver$(OBJ_EXT) libelmersolver$(SHL_EXT)
	$(MPI_F90) $(SH_LINKING_TO_FLAGS) $(RPATH_ELMER) $(B64FLAGS) $(LDFLAGS) -o ElmerSolver_mpi$(EXE_EXT) \
		   Solver$(OBJ_EXT) -L. -lelmersolver $(MPI_LIBS) $(EXTRA_LIBS)
else
ElmerSolver_mpi$(EXE_EXT): Solver$(OBJ_EXT) libelmersolver$(SHL_EXT)
	$(FC) $(SH_LINKING_TO_FLAGS) $(RPATH_ELMER) $(B64FLAGS) $(LDFLAGS) \
	      -o ElmerSolver_mpi$(EXE_EXT) Solver$(OBJ_EXT) \
	      -L. -lelmersolver $(MPI_LIBS) $(EXTRA_LIBS)
endif
endif



ElmerSolver$(EXE_EXT): Solver$(OBJ_EXT) libelmersolver$(SHL_EXT) mpi_stubs$(OBJ_EXT)
	$(FC) $(SH_LINKING_TO_FLAGS) $(RPATH_ELMER) $(B64FLAGS) $(LDFLAGS) \
	      -o ElmerSolver$(EXE_EXT) Solver$(OBJ_EXT) mpi_stubs$(OBJ_EXT) \
	      -L. -lelmersolver $(EXTRA_LIBS) 

ViewFactors$(EXE_EXT): libelmersolver$(SHL_EXT) ViewFactors$(OBJ_EXT) mpi_stubs$(OBJ_EXT)
if USE_WINDOWS_COMPILER
	$(FC) $(SH_LINKING_TO_FLAGS) $(RPATH_ELMER) $(B64FLAGS) $(LDFLAGS) \
	      -o ViewFactors$(EXE_EXT) ViewFactors$(OBJ_EXT) mpi_stubs$(OBJ_EXT) \
	      -L. -lelmersolver viewaxis/libviewaxis.a view3d/libview3d.a $(EXTRA_LIBS) $(BLAS_LIBS)
else
	$(FC) $(SH_LINKING_TO_FLAGS) $(RPATH_ELMER) $(B64FLAGS) $(LDFLAGS) \
	      -o ViewFactors$(EXE_EXT) ViewFactors$(OBJ_EXT) mpi_stubs$(OBJ_EXT) \
	      -L. -lelmersolver viewaxis/libviewaxis.a view3d/libview3d.a $(EXTRA_LIBS) $(BLAS_LIBS) $(STDCXX_LIBS)
endif

GebhardtFactors$(EXE_EXT): libelmersolver$(SHL_EXT) GebhardtFactors$(OBJ_EXT) mpi_stubs$(OBJ_EXT)
	$(FC) $(SH_LINKING_TO_FLAGS) $(RPATH_ELMER) $(B64FLAGS) $(LDFLAGS) \
	      -o GebhardtFactors$(EXE_EXT) GebhardtFactors$(OBJ_EXT) mpi_stubs$(OBJ_EXT) \
	      -L. -lelmersolver viewaxis/libviewaxis.a view3d/libview3d.a $(EXTRA_LIBS)



if USE_WINDOWS_COMPILER
else
ResultToResult$(EXE_EXT): libelmersolver$(SHL_EXT) ResultToResult$(OBJ_EXT)
	$(FC) $(SH_LINKING_TO_FLAGS) $(RPATH_ELMER) $(B64FLAGS) $(LDFLAGS) \
	      -o ResultToResult$(EXE_EXT) ResultToResult$(OBJ_EXT) \
	      -L. -lelmersolver mpi_stubs$(OBJ_EXT) $(EXTRA_LIBS)
endif

ELMER_SOLVER_DATADIR=$(prefix)/share/elmersolver
install-data-local:
	@echo "Installing solvers to ${ELMER_SOLVER_DATADIR}/lib"

	mkdir -p $(DESTDIR)$(ELMER_SOLVER_DATADIR)/lib
	$(CP) elements.def $(DESTDIR)$(ELMER_SOLVER_DATADIR)/lib
	$(CP) SOLVER.KEYWORDS $(DESTDIR)$(ELMER_SOLVER_DATADIR)/lib

	@echo "Installing *.mod to ${ELMER_SOLVER_DATADIR}/include"
	mkdir -p $(DESTDIR)$(ELMER_SOLVER_DATADIR)/include
	$(CP) *.mod $(DESTDIR)$(ELMER_SOLVER_DATADIR)/include


install-exec-local:
	@echo "Installing stuff to ${prefix}/bin"
	$(MKDIR) $(DESTDIR)$(prefix)/bin
	$(CP) ElmerSolver$(EXEEXT) $(DESTDIR)$(prefix)/bin
	$(CP) ViewFactors$(EXEEXT) $(DESTDIR)$(prefix)/bin
	$(CP) GebhardtFactors$(EXEEXT) $(DESTDIR)$(prefix)/bin
if USE_WINDOWS_COMPILER
	$(CP) libelmersolver*$(SHL_EXT) $(DESTDIR)$(prefix)/lib
else
	$(CP) libelmersolver-$(VERSION)$(SHL_EXT) $(DESTDIR)$(prefix)/lib
	$(LN) -sf libelmersolver-$(VERSION)$(SHL_EXT)\
              $(DESTDIR)$(prefix)/lib/libelmersolver$(SHL_EXT)
endif
	$(CP) elmerf90 elmerf90-nosh elmerld $(DESTDIR)$(prefix)/bin
if USE_MPI
	$(CP) ElmerSolver_mpi$(EXEEXT) $(DESTDIR)$(prefix)/bin
endif

clean-local:
	$(RM) $(ALL) *.o *.mod *.f90 *.F *.fpp *~ *.a *.so *.dll \#* mpif mpif.h ElmerSolver ElmerSolver_mpi GebhardtFactors ViewFactors


bin_PROGRAMS = SC2Elmer 
SC2Elmer_SOURCES = SC2Elmer.c


if USE_LOCAL_MPIF_H
mpif: mpif_stub.h
	$(CP) mpif_stub.h mpif.h
	$(CP) mpif_stub.h mpif
else 
mpif: mpif_stub.h
	$(CP) mpif_stub.h mpif
endif

EXTRA_DIST = \
	Adaptive.src \
	AdvectionDiffusion.src \
	AdvectionReaction.src \
	ArteryOutlet.src \
	ArtificialCompressibility.src \
	BandMatrix.src \
	BandwidthOptimize.src \
	BlockSolve.src \
        cholmod.c	      \
	ClusteringMethods.src \
	CoordinateSystems.src \
	CPUTime.c \
	CRSMatrix.src \
	DataToFieldSolver.src \
	DCRComplexSolve.src \
	DefUtils.src \
	DFTSolver.src \
	Differentials.src \
	DiffuseConvectiveAnisotropic.src \
	DiffuseConvectiveGeneralAnisotropic.src \
	DirectSolve.src \
	DistanceSolve.src \
	DivergenceSolver.src \
	EigenSolve.src \
	ElasticSolve.src \
	ElectricForce.src \
	Electrokinetics.src \
	ElementDescription.src \
	ElementUtils.src \
	EliminateDirichlet.src \
	EliminatePeriodic.src \
	ElmerSolver.src \
	EnergyRelease.src \
	FDiffusion3D.src \
	FilterTimeSeries.src \
	ExchangeCorrelations.src \
	fft.c \
	FDiffusion.src \
	Feti.src \
	FindOptimum.src \
	FlowSolve.src \
	FluidicForce.src \
	FluxSolver.src \
	FreeSurface.src \
	FreeSurfaceReduced.src \
	FreeSurfaceSolver.src \
	PhaseChangeSolve.src \
	GebhardtFactors.src \
	GeneralUtils.src \
	HashTable.src \
	HeatSolve.src \
	HelmholtzBEM.src \
	HelmholtzSolve.src \
	Integration.src \
	Interpolation.src \
	Isosurface.src \
	IterSolve.src \
	IterativeMethods.src \
	KESolver.src \
	LevelSet.src \
	LinearAlgebra.src \
	ListMatrix.src \
	Lists.src \
	Load.c \
	LUDecomposition.src \
	MagneticSolve.src \
	MagneticW1Solve.src \
	MainUtils.src \
	MaterialModels.src \
	MaxwellAxiS.src \
	Maxwell.src \
	MaxwellGeneral.src \
	MeshSolve.src \
	MeshUtils.src \
	Messages.src \
	MagnetoDynamics.src \
	MGPrec.src \
	ModelDescription.src \
	MovingElstatSolver.src \
	mpif_stub.h \
	mpi_stubs.c \
	Multigrid.src \
	NavierStokesCylindrical.src \
	NavierStokes.src \
	NavierStokesGeneral.src \
	NonphysicalMeshSolve.src \
	NormalSolver.src \
	ParallelEigenSolve.src \
	ParticleAdvector.src \
	ParticleDynamics.src \
	ParticleUtils.src \
	PElementBase.src \
	PElementMaps.src \
	PoissonBEM.src \
	PoissonBoltzmannSolve.src \
	PhaseChangeSolve.src \
	Radiation.src \
	RadiationFactors.src \
	RateOfChange.src \
	ReloadData.src \
	ReloadInput.src \
	ResultOutputSolve.src \
	ResultToPost.src \
	ResultToResult.src \
	ReynoldsSolver.src \
	RichardsSolver.src \
	RigidBodyReduction.src \
	RigidMeshMapper.src \
	SaveData.src \
	SaveGridData.src \
	ScalarPotentialSolver.src \
	ShearrateSolver.src \
	ShellSolve.src \
	Smitc.src \
	Smoothers.src \
	SolveBand.src \
	solve_cmplx.src \
	solve_real.src \
	Solver.src \
	SolverUtils.src \
	SolveSBand.src \
	SolveHypre.c   \
	SolveTrilinos.cxx   \
	SolveSuperLU.c \
	Spalart-Allmaras.src\
	SParIterComm.src \
	SParIterGlobals.src \
	SParIterPrecond.src \
	SParIterSolver.src \
	SSTKomega.src \
	StatCurrentSolve.src \
	StatElecBoundary.src \
	StatElecSolve.src \
	StatMagSolve.src \
	SteadyPhaseChange.src \
	Stokes.src \
	StreamSolver.src \
	Stress.src \
	StressGeneral.src \
	StressSolve.src \
	StucturedFlowLine.src \
	StructuredMeshMapper.src \
	ThermoElectricSolver.src \
	TimeIntegrate.src \
	TransportEquation.src \
	TransientPhaseChange.src \
	Types.src \
	V2FSolver.src \
	Komega.src \
	SSTKomega.src \
	VankaCreate.src \
	ViewFactors.src \
	VorticitySolver.src \
	Walls.src \
	iso_varying_string.src \
	elements.def \
	f_stubs.c    \
	elmerf90.in  \
	elmerld.in  \
	elmerf90-nosh.in  \
	umf4_f77wrapper.c \
	SOLVER.KEYWORDS
