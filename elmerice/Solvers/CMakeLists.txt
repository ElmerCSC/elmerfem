#SET(WITH_GridDataReader FALSE CACHE BOOL "Include ElmerIce GridDataReader")
SET(WITH_ScatteredDataInterpolator FALSE CACHE BOOL "Include ElmerIce ScatteredDataInterpolator")
#MARK_AS_ADVANCED(WITH_GridDataReader)
MARK_AS_ADVANCED(WITH_ScatteredDataInterpolator)

# ---------------------- #
# -- NetCDF libraries -- #
#----------------------- #
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "Elmer/Ice package: Looking for [NetCDF] & [NetCDF Fortran] libraries")

FIND_PACKAGE(NetCDF)

IF(NetCDF_FOUND)

  SET(HAVE_NETCDF TRUE)
  MARK_AS_ADVANCED(HAVE_NETCDF)
  INCLUDE_DIRECTORIES(${NetCDF_INCLUDE_DIR})
  ADD_DEFINITIONS(-DHAVE_NETCDF)

  MESSAGE(STATUS "  NetCDF:           " "${NetCDF_FOUND}")
  MESSAGE(STATUS "  NetCDF_INC:       " "${NetCDF_INCLUDE_DIR}")
  MESSAGE(STATUS "  NetCDF_LIBS:     " "${NetCDF_LIBRARIES}")
ELSE()
  MESSAGE(STATUS "Library not found: >NetCDF_FOUND< ")
  MESSAGE(WARNING "  \n Missing: >NetCDF_INCLUDE_DIR< , >NetCDF_LIBRARY<,  >NetCDFF_LIBRARY< \n some functionalities will be disabled")
ENDIF(NetCDF_FOUND)

# ---------------------- #
# -- HDF5 libraries -- #
#----------------------- #
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "Elmer/Ice package: Looking for [HDF5] & [HDF5 Fortran] libraries")

FIND_PACKAGE(PHDF5)

IF(PHDF5_FOUND)

  SET(HAVE_HDF5 TRUE)
  MARK_AS_ADVANCED(HAVE_HDF5)
  INCLUDE_DIRECTORIES(${PHDF5_INCLUDE_DIR})
  ADD_DEFINITIONS(-DHAVE_HDF5)

  MESSAGE(STATUS "  HDF5:           " "${PHDF5_FOUND}")
  MESSAGE(STATUS "  HDF5_INC:       " "${PHDF5_INCLUDE_DIR}")
  MESSAGE(STATUS "  HDF5_LIBS:     " "${PHDF5_LIBRARIES}")
ELSE()
  MESSAGE(STATUS "Library not found: >HDF5_FOUND< ")
  MESSAGE(WARNING "  \n Missing: >HDF5_INCLUDE_DIR< , >HDF5_LIBRARY<,  >HDF5F_LIBRARY< \n some functionalities will be disabled")
ENDIF(PHDF5_FOUND)


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
    ./AdjointSSA/AdjointSSA_AdjointSolver.F90 ./AdjointSSA/AdjointSSA_CostDiscSolver.F90
    ./AdjointSSA/AdjointSSA_CostRegSolver.F90 ./AdjointSSA/AdjointSSA_SSASolver.F90
    ./AdjointSSA/AdjointSSA_CostContSolver.F90 ./AdjointSSA/AdjointSSA_CostFluxDivSolver.F90
    ./AdjointSSA/AdjointSSA_GradientSolver.F90 Calving.F90 FrontDisplacement.F90 
    TwoMeshes.F90 ProjectCalving.F90 ComputeCalvingNormal.F90 
    CalvingGeometry.F90 Calving3D.F90 CalvingRemesh.F90
    CalvingFrontAdvance3D.F90 Emergence.F90 SSAmask.F90 
    GlaDSCoupledSolver.F90 GlaDSchannelSolver.F90 Flotation.F90 )

IF(HAVE_MMG)
  SET(ElmerIce_SRC ${ElmerIce_SRC} CalvingRemeshMMG.F90)
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
TARGET_LINK_LIBRARIES(ElmerIceSolvers elmersolver)

IF(HAVE_NETCDF)
TARGET_LINK_LIBRARIES(ElmerIceSolvers ${NetCDF_LIBRARIES})
ENDIF()
IF(HAVE_HDF5)
TARGET_LINK_LIBRARIES(ElmerIceSolvers ${PHDF5_LIBRARIES})
ENDIF()
# IF(HAVE_MMG)
#   TARGET_LINK_LIBRARIES(ElmerIceSolvers ${MMG_LIBRARY})
#   TARGET_LINK_LIBRARIES(ElmerIceSolvers ElmerIce_MeshAdapt3D)
# ENDIF()

INSTALL(TARGETS ElmerIceSolvers LIBRARY DESTINATION "share/elmersolver/lib"
	    RUNTIME DESTINATION "share/elmersolver/lib")

IF(HAVE_NETCDF)
  ADD_SUBDIRECTORY(GridDataReader)
ENDIF()
IF(WITH_ScatteredDataInterpolator)	
  ADD_SUBDIRECTORY(ScatteredDataInterpolator)
ENDIF()
