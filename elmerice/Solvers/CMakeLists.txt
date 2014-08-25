SET(CMAKE_Fortran_MODULE_DIRECTORY
  ${PROJECT_BINARY_DIR}/fmodules CACHE PATH "Directory for Fortran modules")

SET(ElmerIce_SRC AIFlowSolve_nlD2.F90 AIFlowSolve_nlS2.F90 
    CaffeSolver.F90 ComputeDevStressNS.F90 ComputeEigenValues.F90 
    ComputeNormal.F90 ComputeStrainRate.F90 DeformationalHeat.F90 
    EPLSolver.F90 fAandfB_in.F90 FabricSolve.F90 Flowdepth.F90 
    ForceToStress.F90 GetHydrostaticLoads.F90  GolfLaw.F90 
    GroundedSolver.F90  IntegratedVelocity.F90  IDSSolver.F90 
    PorousSolve.F90  pointwise.F90  SIASolver.F90  SSASolver.F90  
    ThicknessSolver.F90  TemperateIce.F90  ExportVertically.F90  
    AdjointSolver.F90  DJDBeta_Adjoint.F90  DJDmu_Adjoint.F90  
    CostSolver_Adjoint.F90  DJDBeta_Robin.F90  DJDmu_Robin.F90  
    CostSolver_Robin.F90  m1qn3.F  Grid2DInterpolator.F90  
    Optimize_m1qn3Parallel.F90  OutputStrainHeating.F90  UpdateExport.F90
    IntegrateVertically.F90 EnthalpySolver.F90)

ADD_LIBRARY(ElmerIceSolvers SHARED ${ElmerIce_SRC})

# Library object
GET_TARGET_PROPERTY(CURRENT_LINK_FLAGS ElmerIceSolvers LINK_FLAGS) 
GET_TARGET_PROPERTY(CURRENT_COMPILE_FLAGS ElmerIceSolvers COMPILE_FLAGS)
IF(CURRENT_LINK_FLAGS)
  SET(CURR_LFLAGS "${CURRENT_LINK_FLAGS};${MPI_Fortran_LINK_FLAGS}")
ELSE()
  SET(CURR_LFLAGS "${MPI_Fortran_LINK_FLAGS}")
ENDIF()
IF(CURRENT_COMPILE_FLAGS)
  SET(CURR_CFLAGS "${CURRENT_COMPILE_FLAGS};${MPI_Fortran_COMPILE_FLAGS}")
ELSE()
  SET(CURR_CFLAGS "${MPI_Fortran_COMPILE_FLAGS}")
ENDIF()
  
SET_TARGET_PROPERTIES(ElmerIceSolvers PROPERTIES LINK_FLAGS "${CURR_LFLAGS}")
SET_TARGET_PROPERTIES(ElmerIceSolvers PROPERTIES COMPILE_FLAGS "${CURR_CFLAGS}")
TARGET_LINK_LIBRARIES(ElmerIceSolvers elmersolver)

INSTALL(TARGETS ElmerIceSolvers LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}")
