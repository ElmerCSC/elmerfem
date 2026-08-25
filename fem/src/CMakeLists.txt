# For the installer process we need a single directory
# where we put all module files
SET(CMAKE_Fortran_MODULE_DIRECTORY
  ${PROJECT_BINARY_DIR}/fmodules CACHE PATH "Directory for Fortran modules")

ADD_SUBDIRECTORY(binio)
ADD_SUBDIRECTORY(modules)
ADD_SUBDIRECTORY(view3d)
ADD_SUBDIRECTORY(viewaxis)

IF(WITH_MPI)
  SET(ELMERSOLVER_OUTPUT_NAME "ElmerSolver_mpi")
ELSE()
  SET(ELMERSOLVER_OUTPUT_NAME "ElmerSolver")
ENDIF()

SET(prefix "${CMAKE_INSTALL_PREFIX}")

SET(solverlib_SOURCES AddrFunc.F90 NavierStokes.F90 NavierStokesGeneral.F90
  NavierStokesCylindrical.F90 Lists.F90
  DiffuseConvectiveAnisotropic.F90 LoadMod.F90
  DiffuseConvectiveGeneralAnisotropic.F90 PElementMaps.F90
  PElementBase.F90 ElementBasis.F90 ElementGeometry.F90 ElemInfo.F90 ElementDescription.F90 Integration.F90 ListMatrixArray.F90
  ModelDescription.F90 GeneralUtils.F90 Constitutive.F90 Stress.F90 StressGeneral.F90 ModelLumping.F90
  LinearAlgebra.F90 CoordinateSystems.F90 ListMatrix.F90 CRSMatrix.F90
  BandMatrix.F90 BandwidthOptimize.F90 BlockSolve.F90
  MaterialModels.F90 DirectSolve.F90 IterSolve.F90
  IterativeMethods.F90 TimeIntegrate.F90 Types.F90 SolveBand.F90
  ElementUtils.F90 Radiation.F90 fft.c Load.c Differentials.F90
  FreeSurface.F90
  Walls.F90 MatrixScaling.F90 SolverBasics.F90 BoundaryConditionUtils.F90 SolveCore.F90 ContactUtils.F90 ProjectorUtils.F90 SolverUtils.F90 SolveSBand.F90 CPUTime.c Interpolation.F90
  MainUtils.F90 Adaptive.F90 EigenSolve.F90 HashTable.F90 MortarUtils.F90
  MeshAllocations.F90 MeshIO.F90 GeometryFitting.F90
  MeshBasics.F90 MeshGraph.F90 MeshTransform.F90 MeshSplit.F90 MeshTagging.F90 MeshLoad.F90 MeshExtrusion.F90 MeshUtils.F90 SaveUtils.F90 SParIterGlobals.F90 SParIterComm.F90
  SParIterPrecond.F90 SParIterSolver.F90 Messages.F90 Multigrid.F90
  Smoothers.F90 ClusteringMethods.F90 ParallelUtils.F90
  ParallelEigenSolve.F90 solve_cmplx.F90 solve_real.F90 MGPrec.F90
  DefUtils.F90 MeshPartition.F90 MeshRemeshing.F90 CutFEMUtils.F90
  RadiationFactors.F90 f_stubs.c GFortranLocaleFix.c
  ExchangeCorrelations.F90 SolveHypre.c amgx.c rocalution.cpp SolverActivate_x.F90
  SolveSuperLU.c
  umf4_f77wrapper.c VankaCreate.F90 ParticleUtils.F90 Feti.F90
  cholmod.c InterpolateMeshToMesh.F90 InterpVarToVar.F90 PElementVisual.F90 IpFieldInterface.F90
  LinearForms.F90 H1Basis.F90 CircuitUtils.F90 BackwardError.F90
  ElmerSolver.F90 MagnetoDynamicsUtils.F90 VectorHelmholtzUtils.F90 LumpingUtils.F90
  ZirkaHysteresis.F90 SolidMechanicsUtils.F90 MatrixAssembly.F90
  OptimizationUtils.F90 MeshGenerate.F90 ViewUtils.F90)

#SET_PROPERTY(SOURCE MaxwellAxiS.F90 PROPERTY
#COMPILE_DEFINITIONS FULL_INDUCTION)

IF(WITH_EXTOPTIM)
  LIST(APPEND solverlib_SOURCES minpack.F90 minpack_kinds.F90
  PowellOptBobyqa.F90 PowellOptNewuoa.F90)
  SET(HAVE_EXTOPTIM TRUE)
ENDIF()

IF(WITH_LUA)
  INCLUDE_DIRECTORIES(${LUA_INCLUDE_DIR})
  LIST(APPEND solverlib_SOURCES Lua.F90 elmer_lua_iface.c)
  INSTALL(DIRECTORY lua-scripts DESTINATION share/elmersolver/
    PATTERN lua-scripts/*.lua)
  FILE(GLOB DEFAULT_LUA_SCRIPTS "lua-scripts/*.lua")
  FILE(COPY ${DEFAULT_LUA_SCRIPTS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/share/elmersolver/lua-scripts/)
ELSE()
  LIST(APPEND solverlib_SOURCES Lua.F90)
ENDIF()

IF(HAVE_FETI4I)
  LIST(APPEND solverlib_SOURCES ${FETI4I_INTERFACE_SOURCE})
ENDIF()

IF(HAVE_ADIOS2)
  LIST(APPEND solverlib_SOURCES ADIOS2Utils.F90)
ENDIF(HAVE_ADIOS2)

#FILE(GLOB SRC_FILES *.src)
#FOREACH(FNAME ${SRC_FILES})
#GET_FILENAME_COMPONENT(BASENAME ${FNAME} NAME_WE)
#ADD_CUSTOM_COMMAND(
  #OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.F90
  #COMMAND ${CMAKE_COMMAND} -E copy ${FNAME}
    #${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.F90
  #DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${BASENAME}.src
  #)
#ENDFOREACH()

INCLUDE_DIRECTORIES("${CMAKE_BINARY_DIR}/fem/src")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/fhutiter/src")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}/binio")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/fhutiter/src")

# Extract additional compile flags
GET_DIRECTORY_PROPERTY( ELMER_CDEFS COMPILE_DEFINITIONS )
FOREACH( d ${ELMER_CDEFS} )
 SET(ELMER_F90FLAGS "${ELMER_F90FLAGS} -D${d}")
ENDFOREACH()


# Extract compile flags arising from build type
IF(CMAKE_BUILD_TYPE)
  STRING(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPCASE)
  SET(ELMER_F90FLAGS "${ELMER_F90FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
  MARK_AS_ADVANCED(CMAKE_BUILD_TYPE_UPCASE)
ENDIF(CMAKE_BUILD_TYPE)

# Build elmerf90 and elmerld as compiled C binaries (replaces the shell script / .bat)
IF(WITH_ElmerIce)
  SET(HAVE_ELMERICE_INT 1)
  set(ELMERF90_ELMERICESOLVERS_LIB "\$\<TARGET_FILE_NAME:ElmerIceSolvers>")
  set(ELMERF90_ELMERICEUSF_LIB "\$\<TARGET_FILE_NAME:ElmerIceUSF>")
  set(ELMERF90_ELMERICEUTILS_LIB "\$\<TARGET_FILE_NAME:ElmerIceUtils>")
ELSE()
  SET(HAVE_ELMERICE_INT 0)
  set(ELMERF90_ELMERICESOLVERS_LIB "")
  set(ELMERF90_ELMERICEUSF_LIB "")
  set(ELMERF90_ELMERICEUTILS_LIB "")
ENDIF()
IF(HAVE_MMG)
  SET(HAVE_MMG_INT 1)
ELSE()
  SET(HAVE_MMG_INT 0)
ENDIF()
IF(HAVE_PARMMG)
  SET(HAVE_PARMMG_INT 1)
ELSE()
  SET(HAVE_PARMMG_INT 0)
ENDIF()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/elmerf90-helper.in.h
  ${CMAKE_CURRENT_BINARY_DIR}/elmerf90-helper.out.h @ONLY)
FILE(GENERATE
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/elmerf90-helper.h
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/elmerf90-helper.out.h)

ADD_EXECUTABLE(elmerf90_bin ${CMAKE_CURRENT_SOURCE_DIR}/elmerf90.c)
SET_TARGET_PROPERTIES(elmerf90_bin PROPERTIES OUTPUT_NAME "elmerf90")

ADD_EXECUTABLE(elmerld_bin ${CMAKE_CURRENT_SOURCE_DIR}/elmerld.c)
SET_TARGET_PROPERTIES(elmerld_bin PROPERTIES OUTPUT_NAME "elmerld")

# Create pkg-config file if requested
IF(CREATE_PKGCONFIG_FILE)
  CONFIGURE_FILE(elmer.in.pc elmer.pc @ONLY)
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/elmer.pc
    DESTINATION "share/pkgconfig")
ENDIF(CREATE_PKGCONFIG_FILE)

# Copy elements.def and SOLVER.KEYWORDS (to enable testing)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/elements.def
  ${CMAKE_CURRENT_BINARY_DIR}/elements.def COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/SOLVER.KEYWORDS
  ${CMAKE_CURRENT_BINARY_DIR}/SOLVER.KEYWORDS COPYONLY)

IF (NOT MPI_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mpif_stub.h ${CMAKE_CURRENT_BINARY_DIR}/mpif.h COPYONLY)
  # add path to stub header to dummy import target
  target_include_directories(Elmer::MPI_Fortran INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()

ADD_LIBRARY(mpi_stubs SHARED mpif_stubs.F90)
ADD_LIBRARY(elmersolver SHARED ${solverlib_SOURCES})
SET_TARGET_PROPERTIES(mpi_stubs
                      PROPERTIES LINKER_LANGUAGE Fortran)
SET_TARGET_PROPERTIES(elmersolver
                      PROPERTIES LINKER_LANGUAGE Fortran)

#We want MMG & Zoltan in the RPATH (if found), but given how many other
#components use ELMERSOLVER_RPATH_STRING, it seemed best to make
#a local copy just for ElmerSolver and libelmersolver
SET(ELMERSOLVER_RPATH_STRING_MOD "${ELMERSOLVER_RPATH_STRING}")
IF(HAVE_MMG)
  TARGET_LINK_LIBRARIES(elmersolver ${MMG_LIBRARY})
  SET(ELMERSOLVER_RPATH_STRING_MOD "${ELMERSOLVER_RPATH_STRING_MOD}:${MMG_LIBDIR}")
  SET(ELMERLIB_RPATH_STRING "${ELMERLIB_RPATH_STRING}:${MMG_LIBDIR}")
  #TODO - add to RPath
ENDIF()

IF(HAVE_PARMMG)
  TARGET_LINK_LIBRARIES(elmersolver ${PARMMG_LIBRARY})
  SET(ELMERSOLVER_RPATH_STRING_MOD "${ELMERSOLVER_RPATH_STRING_MOD}:${PARMMG_LIBDIR}")
  SET(ELMERLIB_RPATH_STRING "${ELMERLIB_RPATH_STRING}:${PARMMG_LIBDIR}")
  #TODO - add to RPath
ENDIF()

IF(HAVE_ZOLTAN)
  IF(USE_SYSTEM_ZOLTAN)
    TARGET_LINK_LIBRARIES(elmersolver Zoltan::zoltan)
  ELSE()
    ADD_DEPENDENCIES(elmersolver zoltan)
    TARGET_LINK_LIBRARIES(elmersolver ${ZOLTAN_LIBRARY})
  ENDIF()
  IF(Zoltan_ENABLE_ParMETIS)
    TARGET_LINK_LIBRARIES(elmersolver ${ParMetis_LIBRARIES})
  ENDIF()
  # SET(ELMERSOLVER_RPATH_STRING_MOD "${ELMERSOLVER_RPATH_STRING_MOD}/:${ZOLTAN_LIBDIR}")
  # SET(ELMERLIB_RPATH_STRING "${ELMERLIB_RPATH_STRING}/:${ZOLTAN_LIBDIR}")
ENDIF()

# ElmerSolver libraries
SET(ELMERSOLVER_LIBRARIES matc fhuti binio arpack)

IF(UMFPACK_FOUND)
  # external implementation of UMFPACK
  LIST(APPEND ELMERSOLVER_LIBRARIES ${UMFPACK_LIBRARIES})

  TARGET_INCLUDE_DIRECTORIES(elmersolver PRIVATE ${UMFPACK_INCLUDE_DIR})
  #LIST(APPEND ELMERSOLVER_LIBRARIES ${UMFPACK_LIBRARIES})
  #TARGET_LINK_LIBRARIES(elmersolver "${UMFPACK_LIBRARIES}")
ELSE()
  # bundled UMFPACK
  LIST(APPEND ELMERSOLVER_LIBRARIES umfpack amd)
  TARGET_INCLUDE_DIRECTORIES(elmersolver PRIVATE ${PROJECT_SOURCE_DIR}/umfpack/src/umfpack/include)
ENDIF()

LIST(APPEND ELMERSOLVER_LIBRARIES
  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
  ${CMAKE_DL_LIBS})

# MUMPS: add includes and libraries regardless of MPI (supports serial build)
IF(Mumps_LIBRARIES)
  INCLUDE_DIRECTORIES(${Mumps_INCLUDE_DIR})
  LIST(APPEND ELMERSOLVER_LIBRARIES ${Mumps_LIBRARIES})
ENDIF()

IF(WITH_LUA)
  LIST(APPEND ELMERSOLVER_LIBRARIES ${LUA_LIBRARIES})
ENDIF()

IF(HAVE_ADIOS2)
  IF(WITH_MPI)
    TARGET_LINK_LIBRARIES(elmersolver adios2::fortran_mpi adios2::fortran)
  ELSE()
    TARGET_LINK_LIBRARIES(elmersolver adios2::fortran)
  ENDIF()
ENDIF(HAVE_ADIOS2)

ADD_EXECUTABLE(Solver_TGT Solver.F90)
SET_TARGET_PROPERTIES(Solver_TGT PROPERTIES OUTPUT_NAME "${ELMERSOLVER_OUTPUT_NAME}")

IF(WIN32)
  # Increase stack size from the default 1 MB to 16 MB.
  TARGET_LINK_OPTIONS(Solver_TGT PRIVATE
    $<IF:$<Fortran_COMPILER_ID:MSVC>,/STACK:0x1000000,-Wl$<COMMA>--stack$<COMMA>0x1000000>)
ELSE()
  MESSAGE(STATUS " ELMERSOLVER_RPATH_STRING_MOD " ${ELMERSOLVER_RPATH_STRING_MOD})
  MESSAGE(STATUS " ELMERLIB_RPATH_STRING " ${ELMERLIB_RPATH_STRING})

   SET_TARGET_PROPERTIES(Solver_TGT PROPERTIES INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING_MOD}")
   SET_TARGET_PROPERTIES(elmersolver PROPERTIES INSTALL_RPATH "${ELMERLIB_RPATH_STRING}")
ENDIF()

IF(NOT(WITH_MPI))
  TARGET_LINK_LIBRARIES(elmersolver mpi_stubs ${ELMERSOLVER_LIBRARIES})
  TARGET_LINK_LIBRARIES(Solver_TGT elmersolver)
ENDIF()

IF(WITH_MPI)
  #ADD_LIBRARY(elmersolver_mpi SHARED ${solverlib_SOURCES})

  IF(MKL_CPARDISO_FOUND AND NOT Mumps_LIBRARIES)
    # SCALAPACK is needed with CPARDISO when MUMPS is not present
    LIST(APPEND ELMERSOLVER_LIBRARIES ${SCALAPACK_LIBRARIES})
  ENDIF()

  IF(Hypre_LIBRARIES)
    LIST(APPEND ELMERSOLVER_LIBRARIES Elmer::MPI_C ${Hypre_LIBRARIES})
  ENDIF()

  IF(XIOS_LIBRARIES)
    # XIOS depends on NetCDF (mandatory) and HDF5 (optional), since linking order matters, XIOS must be added before the other two
    LIST(APPEND ELMERSOLVER_LIBRARIES ${XIOS_LIBRARIES})
  ENDIF()

  IF(NETCDF_LIBRARIES)
    LIST(APPEND ELMERSOLVER_LIBRARIES ${NETCDF_LIBRARIES})
    IF(PHDF5_LIBRARIES)
      LIST(APPEND ELMERSOLVER_LIBRARIES ${PHDF5_LIBRARIES})
    ENDIF()
  ENDIF()

  IF(FETI4I_LIBRARIES)
    LIST(APPEND ELMERSOLVER_LIBRARIES ${FETI4I_LIBRARIES})
  ENDIF()

  # Executable
  TARGET_LINK_LIBRARIES(Solver_TGT Elmer::MPI_Fortran elmersolver ${ELMERSOLVER_LIBRARIES})

  # Library object
  TARGET_LINK_LIBRARIES(elmersolver Elmer::MPI_Fortran ${ELMERSOLVER_LIBRARIES} parpack)

  SET(SERIAL_SOLVER_FILEPATH "${CMAKE_INSTALL_PREFIX}/bin/ElmerSolver")
  SET(MPI_SOLVER_FILEPATH "${CMAKE_INSTALL_PREFIX}/bin/${ELMERSOLVER_OUTPUT_NAME}")
ENDIF()

IF(MPI_Fortran_HAVE_F90_MODULE)
  TARGET_COMPILE_DEFINITIONS(elmersolver PRIVATE ELMER_HAVE_MPI_MODULE)
ELSE()
  # Either the header from an actual MPI implementation or the replacement stub
  TARGET_COMPILE_DEFINITIONS(elmersolver PRIVATE ELMER_HAVE_MPIF_HEADER)
ENDIF()
IF(WITH_ROCALUTION)
  TARGET_LINK_LIBRARIES(elmersolver roc::rocalution)
ENDIF()

IF(WITH_AMGX)
  TARGET_LINK_LIBRARIES(elmersolver Elmer::MPI_C "${AMGX_LIBRARY};${CUDA_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(elmersolver PRIVATE "${CUDA_INCLUDE_DIRS}")
ENDIF()

IF (WITH_CHOLMOD)
  SET(CURR_CFLAGS "-DHAVE_CHOLMOD")
  TARGET_INCLUDE_DIRECTORIES(elmersolver PRIVATE "${CHOLMOD_INCLUDE_DIR};${SPQR_INCLUDE_DIR}")
  TARGET_LINK_LIBRARIES(elmersolver "${CHOLMOD_LIBRARIES};${SPQR_LIBRARIES}")
ENDIF()

# elmersolver -library
#TARGET_LINK_LIBRARIES(elmersolver ${ELMERSOLVER_LIBRARIES})

## ElmerSolver serial
#IF(NOT(NO_SERIAL_BINARY))
  #TARGET_LINK_LIBRARIES(elmersolver mpi_stubs ${ELMERSOLVER_LIBRARIES})
  #TARGET_LINK_LIBRARIES(ElmerSolver elmersolver)
#ENDIF()

# ElmerSolver
#ADD_EXECUTABLE(ElmerSolver Solver.F90)
#TARGET_LINK_LIBRARIES(ElmerSolver ${ELMERSOLVER_LIBRARIES}
                                  #mpi_stubs)

# ViewFactors
ADD_EXECUTABLE(ViewFactors ViewFactors.F90)
SET_TARGET_PROPERTIES(ViewFactors PROPERTIES LINKER_LANGUAGE Fortran)
# TARGET_INCLUDE_DIRECTORIES(ViewFactors PUBLIC ${CMAKE_BINARY_DIR}/fem/src/view3d)
# TARGET_INCLUDE_DIRECTORIES(ViewFactors PUBLIC ${CMAKE_BINARY_DIR}/fem/src/viewfactors)
# FIND_LIBRARY(VIEW3D_LIBRARIES PATHS ${CMAKE_CURRENT_BINARY_DIR}/view3d)
# FIND_LIBRARY(VIEWAXIS_LIBRARIES PATHS ${CMAKE_CURRENT_BINARY_DIR}/viewaxis)
TARGET_LINK_LIBRARIES(ViewFactors view3d viewaxis elmersolver)
INSTALL(TARGETS ViewFactors RUNTIME DESTINATION "bin")

IF(NOT(WIN32))
  SET_TARGET_PROPERTIES(ViewFactors PROPERTIES INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING}")
ENDIF()

# Radiators
ADD_EXECUTABLE(Radiators Radiators.F90)
SET_TARGET_PROPERTIES(Radiators PROPERTIES LINKER_LANGUAGE Fortran)
# TARGET_INCLUDE_DIRECTORIES(Radiators PUBLIC ${CMAKE_BINARY_DIR}/fem/src/view3d)
# TARGET_INCLUDE_DIRECTORIES(Radiators PUBLIC ${CMAKE_BINARY_DIR}/fem/src/viewfactors)
# FIND_LIBRARY(VIEW3D_LIBRARIES PATHS ${CMAKE_CURRENT_BINARY_DIR}/view3d)
# FIND_LIBRARY(VIEWAXIS_LIBRARIES PATHS ${CMAKE_CURRENT_BINARY_DIR}/viewaxis)
TARGET_LINK_LIBRARIES(Radiators view3d viewaxis elmersolver)
INSTALL(TARGETS Radiators RUNTIME DESTINATION "bin")

IF(NOT(WIN32))
  SET_TARGET_PROPERTIES(Radiators PROPERTIES INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING}")
ENDIF()


# Installation rules
INSTALL(FILES elements.def SOLVER.KEYWORDS DESTINATION "share/elmersolver/lib")


INSTALL(TARGETS Solver_TGT elmersolver
  RUNTIME DESTINATION "bin"
  LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${ELMER_INSTALL_LIB_DIR}")
IF(NOT WITH_MPI)
  INSTALL(TARGETS mpi_stubs
    RUNTIME DESTINATION "bin"
    LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}"
    ARCHIVE DESTINATION "${ELMER_INSTALL_LIB_DIR}")
ENDIF()

INSTALL(TARGETS elmerf90_bin RUNTIME DESTINATION "bin")
INSTALL(TARGETS elmerld_bin RUNTIME DESTINATION "bin")


INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/.
                  DESTINATION "share/elmersolver/include")

IF(WITH_MPI)
  IF(NOT(WIN32))
    # Unix: expose the generic name "ElmerSolver" as a symlink to the actual
    # MPI binary "ElmerSolver_mpi" (cheap and unprivileged here). Unchanged.
    ADD_CUSTOM_TARGET(elmersolver_link ALL COMMAND ${CMAKE_COMMAND} -E create_symlink
      ElmerSolver_mpi${CMAKE_EXECUTABLE_SUFFIX}
      ElmerSolver${CMAKE_EXECUTABLE_SUFFIX})
    ADD_DEPENDENCIES(elmersolver_link Solver_TGT)
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/ElmerSolver${CMAKE_EXECUTABLE_SUFFIX}
      DESTINATION "bin")
  ELSE()
    # Windows: COPY rather than symlink. Creating a symlink needs
    # SeCreateSymbolicLinkPrivilege (Developer Mode or an elevated shell),
    # which a normal account lacks, so the previous create_symlink step failed
    # the whole build for ordinary users. A copy works unprivileged (and
    # survives being zipped/packaged). copy_if_different keeps it in sync on
    # rebuilds; it runs after Solver_TGT so the source binary already exists.
    ADD_CUSTOM_COMMAND(TARGET Solver_TGT POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:Solver_TGT>
        ${CMAKE_CURRENT_BINARY_DIR}/ElmerSolver${CMAKE_EXECUTABLE_SUFFIX}
      COMMENT "Copying ${ELMERSOLVER_OUTPUT_NAME} to generic name ElmerSolver")
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/ElmerSolver${CMAKE_EXECUTABLE_SUFFIX}
      DESTINATION "bin")
  ENDIF()
ENDIF()

IF(CMAKE_SYSTEM_NAME MATCHES "Windows" AND CPACK_BUNDLE_EXTRA_WINDOWS_DLLS)
  # Collect list of dependent DLLs using bundledlls
  # This assumes that ElmerFEM is built using MSYS2.
  FIND_PROGRAM(MSYS2_BASH NAMES bash REQUIRED)
  ADD_CUSTOM_TARGET(elmersolver_bundle ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bundle/bin
    COMMAND ${MSYS2_BASH} -lc
      "${CMAKE_CURRENT_SOURCE_DIR}/../../contrib/msys2-bundledlls/bundledlls \
      $<TARGET_FILE:elmersolver> ${CMAKE_BINARY_DIR}/bundle/bin"
  )

  ADD_DEPENDENCIES(elmersolver_bundle elmersolver)

  IF(NOT EXTERNAL_ARPACK)
    ADD_CUSTOM_COMMAND(
      TARGET elmersolver_bundle POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/bundle/bin/libarpack.dll
    )
  ENDIF()

  IF(NOT EXTERNAL_PARPACK)
    ADD_CUSTOM_COMMAND(
      TARGET elmersolver_bundle POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/bundle/bin/libparpack.dll
    )
  ENDIF()

  INSTALL(DIRECTORY ${CMAKE_BINARY_DIR}/bundle/bin/
    DESTINATION "bin")
ENDIF()
