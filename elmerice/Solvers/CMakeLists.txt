#SET(WITH_GridDataReader FALSE CACHE BOOL "Include ElmerIce GridDataReader")
SET(WITH_ScatteredDataInterpolator FALSE CACHE BOOL "Include ElmerIce ScatteredDataInterpolator")
#MARK_AS_ADVANCED(WITH_GridDataReader)
MARK_AS_ADVANCED(WITH_ScatteredDataInterpolator)

# ---------------------- #
# -- LOOK FOR MMG     -- #
#----------------------- #
IF(HAVE_MMG)
  ADD_SUBDIRECTORY(MeshAdaptation_2D)
ENDIF()

SET(CMAKE_Fortran_MODULE_DIRECTORY
  ${PROJECT_BINARY_DIR}/fmodules CACHE PATH "Directory for Fortran modules")

SET(ElmerIce_SRC ElmerIceUtils.F90 AIFlowSolve_nlD2.F90 AIFlowSolve_nlS2.F90 
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
    IntegrateVertically.F90 EnthalpySolver.F90 SubShelfMelt.F90
    ./Adjoint/Adjoint_LinearSolver.F90 ./Adjoint/Adjoint_CostDiscSolver.F90
    ./Adjoint/Adjoint_CostContSolver.F90 ./Adjoint/Adjoint_CostRegSolver.F90 
    ./Adjoint/Adjoint_GradientValidation.F90
    ./AdjointStokes/AdjointStokes_GradientMu.F90 
    ./AdjointStokes/AdjointStokes_GradientBetaSolver.F90  
    ./AdjointSSA/AdjointSSA_AdjointSolver.F90 ./AdjointSSA/AdjointSSA_CostDiscSolver.F90
    ./AdjointSSA/AdjointSSA_CostRegSolver.F90 ./AdjointSSA/AdjointSSA_SSASolver.F90
    ./AdjointSSA/AdjointSSA_CostContSolver.F90 ./AdjointSSA/AdjointSSA_CostFluxDivSolver.F90
    ./AdjointSSA/AdjointSSA_CostTaubSolver.F90
    ./AdjointSSA/AdjointSSA_GradientSolver.F90 
    ./AdjointThickness/AdjointThickness_GradientSolver.F90   ./AdjointThickness/AdjointThickness_ThicknessSolver.F90
    ./Permafrost/PermafrostMaterials.F90 ./Permafrost/Permafrost_Utils.F90 ./Permafrost/Permafrost_HTEQ.F90
    ./Permafrost/Permafrost_Darcy.F90 ./Permafrost/Permafrost_solute.F90 ./Permafrost/Permafrost_solid.F90
    SurfaceBoundaryEnthalpy.F90
    Calving.F90 FrontDisplacement.F90 
    TwoMeshes.F90 ProjectCalving.F90 ComputeCalvingNormal.F90 
    CalvingGeometry.F90 Calving3D.F90 Calving3D_lset.F90 CalvingGlacierAdvance3D.F90 CalvingRemesh.F90
    CalvingFrontAdvance3D.F90 Emergence.F90 SSAmask.F90 
    GlaDSCoupledSolver.F90 GlaDSchannelSolver.F90 Flotation.F90
    BasalMelt3D.F90 CalvingHydroInterp.F90 HydroRestart.F90 
    GMValid.F90 Scalar_OUTPUT_Glacier.F90 IcyMaskSolver.F90
    Weertman2Coulomb.F90)

SET(ElmerIce_SRC ${ElmerIce_SRC} ./Covarianceutils/CovarianceUtils.F90 ./Covarianceutils/BackgroundErrorCostSolver.F90 ./Covarianceutils/CovarianceVectorMultiplySolver.F90 ./Covarianceutils/GaussianSimulationSolver.F90)

IF(HAVE_ODEPack)
  SET(ElmerIce_SRC ${ElmerIce_SRC} PlumeSolver.F90 )
ENDIF()

IF(HAVE_MMG)
  SET(ElmerIce_SRC ${ElmerIce_SRC} CalvingRemeshMMG.F90 ) #MMG3DSolver.F90 <- defunct?
ENDIF()

IF(HAVE_XIOS)
  SET(ElmerIce_SRC ${ElmerIce_SRC} OutPutSolvers/XIOSOutputSolver.F90)
ENDIF()

IF(HAVE_NETCDF)
  SET(ElmerIce_SRC ${ElmerIce_SRC} UGridDataReader.F90)
ENDIF()

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

SET_TARGET_PROPERTIES(ElmerIceSolvers PROPERTIES PREFIX "")  
SET_TARGET_PROPERTIES(ElmerIceSolvers PROPERTIES LINK_FLAGS "${CURR_LFLAGS}")
SET_TARGET_PROPERTIES(ElmerIceSolvers PROPERTIES COMPILE_FLAGS "${CURR_CFLAGS}")
SET_TARGET_PROPERTIES(ElmerIceSolvers PROPERTIES LINKER_LANGUAGE Fortran)
IF(NOT(WIN32))
    SET_TARGET_PROPERTIES(ElmerIceSolvers PROPERTIES
                          INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING};${ELMER_SOLVER_HOME}/lib/")
    SET_TARGET_PROPERTIES(ElmerIceSolvers PROPERTIES
                          INSTALL_RPATH_USE_LINK_PATH TRUE)
ENDIF()
TARGET_LINK_LIBRARIES(ElmerIceSolvers elmersolver ElmerIceUtils)

IF(HAVE_NETCDF)
TARGET_LINK_LIBRARIES(ElmerIceSolvers ${NETCDF_LIBRARIES})
ENDIF()
IF(HAVE_HDF5)
TARGET_LINK_LIBRARIES(ElmerIceSolvers ${PHDF5_LIBRARIES})
ENDIF()
IF(HAVE_XIOS)
  INCLUDE_DIRECTORIES(${XIOS_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES(ElmerIceSolvers ${XIOS_LIBRARIES})
ENDIF()
# IF(HAVE_MMG)
#   TARGET_LINK_LIBRARIES(ElmerIceSolvers ${MMG_LIBRARY})
#   TARGET_LINK_LIBRARIES(ElmerIceSolvers ElmerIce_MeshAdapt3D)
# ENDIF()

IF(HAVE_PROJ)
  INCLUDE_DIRECTORIES(${PROJ_INCLUDE_DIR})
  INCLUDE_DIRECTORIES(${FORTRANGIS_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES(ElmerIceSolvers ${PROJ_LIBRARY})
  TARGET_LINK_LIBRARIES(ElmerIceSolvers ${FORTRANGIS_LIBRARY})
ENDIF()

INSTALL(TARGETS ElmerIceSolvers LIBRARY DESTINATION "share/elmersolver/lib"
	    RUNTIME DESTINATION "share/elmersolver/lib")

IF(HAVE_NETCDF)
  ADD_SUBDIRECTORY(GridDataReader)
ENDIF()
IF(WITH_ScatteredDataInterpolator)	
  ADD_SUBDIRECTORY(ScatteredDataInterpolator)
ENDIF()
