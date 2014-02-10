ADD_SUBDIRECTORY(binio)
ADD_SUBDIRECTORY(modules)

SET(prefix "${CMAKE_INSTALL_PREFIX}")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/elmerf90.in.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/elmerf90)

SET(solverlib_SOURCES NavierStokes.F90 NavierStokesGeneral.F90
  NavierStokesCylindrical.F90 Lists.F90
  DiffuseConvectiveAnisotropic.F90
  DiffuseConvectiveGeneralAnisotropic.F90 PElementMaps.F90
  PElementBase.F90 ElementDescription.F90 Integration.F90
  ModelDescription.F90 GeneralUtils.F90 Stress.F90 StressGeneral.F90
  LinearAlgebra.F90 CoordinateSystems.F90 ListMatrix.F90 CRSMatrix.F90
  BandMatrix.F90 BandwidthOptimize.F90 BlockSolve.F90
  MaterialModels.F90 DirectSolve.F90 IterSolve.F90
  IterativeMethods.F90 TimeIntegrate.F90 Types.F90 SolveBand.F90
  ElementUtils.F90 Radiation.F90 fft.c Load.c Differentials.F90
  FreeSurface.F90 Maxwell.F90 MaxwellAxiS.F90 MaxwellGeneral.F90
  Walls.F90 SolverUtils.F90 SolveSBand.F90 CPUTime.c Interpolation.F90
  MainUtils.F90 Adaptive.F90 EigenSolve.F90 HashTable.F90
  MeshUtils.F90 SParIterGlobals.F90 SParIterComm.F90
  SParIterPrecond.F90 SParIterSolver.F90 Messages.F90 Multigrid.F90
  Smoothers.F90 ClusteringMethods.F90 ParallelUtils.F90
  ParallelEigenSolve.F90 solve_cmplx.F90 solve_real.F90 MGPrec.F90
  DefUtils.F90 LUDecomposition.F90 RadiationFactors.F90 f_stubs.c
  ExchangeCorrelations.F90 SolveHypre.c SolveTrilinos.cxx
  SolveSuperLU.c iso_varying_string.F90 umf4_f77wrapper.c
  VankaCreate.F90 ParticleUtils.F90 Feti.F90 cholmod.c
  InterpolateMeshToMesh.F90 BackwardError.F90 ElmerSolver.F90)

FILE(GLOB SRC_FILES *.src)
FOREACH(FNAME ${SRC_FILES})
GET_FILENAME_COMPONENT(BASENAME ${FNAME} NAME_WE)
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.F90
  COMMAND ${CMAKE_COMMAND} -E copy ${FNAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.F90
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${BASENAME}.src
  )
ENDFOREACH()

INCLUDE_DIRECTORIES("${CMAKE_BINARY_DIR}/fem/src")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/hutiter/src")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}/binio")

ADD_LIBRARY(elmersolver SHARED ${solverlib_SOURCES})
TARGET_LINK_LIBRARIES(elmersolver ${MPI_Fortran_LIBRARIES} matc eioc
  eiof umfpack amd huti ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} arpack
  parpack binio)

ADD_EXECUTABLE(ElmerSolver_mpi Solver.F90)

TARGET_LINK_LIBRARIES(ElmerSolver_mpi elmersolver matc huti eioc eiof
  umfpack amd ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} arpack parpack binio
  ${CMAKE_DL_LIBS})

INSTALL(FILES elements.def SOLVER.KEYWORDS DESTINATION "share/elmersolver/lib")
INSTALL(TARGETS ElmerSolver_mpi elmersolver RUNTIME DESTINATION "bin" LIBRARY DESTINATION "lib") 
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/elmerf90 DESTINATION "bin")
