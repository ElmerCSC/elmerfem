# For the installer process we need a single directory where we put all module
# files
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/fmodules
    CACHE PATH "Directory for Fortran modules")

add_subdirectory(binio)
add_subdirectory(modules)
add_subdirectory(view3d)
add_subdirectory(viewaxis)

if(WITH_MPI)
  set(ELMERSOLVER_OUTPUT_NAME "ElmerSolver_mpi")
else()
  set(ELMERSOLVER_OUTPUT_NAME "ElmerSolver")
endif()

set(prefix "${CMAKE_INSTALL_PREFIX}")

set(solverlib_SOURCES
    AddrFunc.F90
    NavierStokes.F90
    NavierStokesGeneral.F90
    NavierStokesCylindrical.F90
    Lists.F90
    DiffuseConvectiveAnisotropic.F90
    LoadMod.F90
    DiffuseConvectiveGeneralAnisotropic.F90
    PElementMaps.F90
    PElementBase.F90
    ElementDescription.F90
    Integration.F90
    ListMatrixArray.F90
    ModelDescription.F90
    GeneralUtils.F90
    Stress.F90
    StressGeneral.F90
    LinearAlgebra.F90
    CoordinateSystems.F90
    ListMatrix.F90
    CRSMatrix.F90
    BandMatrix.F90
    BandwidthOptimize.F90
    BlockSolve.F90
    MaterialModels.F90
    DirectSolve.F90
    IterSolve.F90
    IterativeMethods.F90
    TimeIntegrate.F90
    Types.F90
    SolveBand.F90
    ElementUtils.F90
    Radiation.F90
    fft.c
    Load.c
    Differentials.F90
    FreeSurface.F90
    Maxwell.F90
    MaxwellAxiS.F90
    MaxwellGeneral.F90
    Walls.F90
    SolverUtils.F90
    SolveSBand.F90
    CPUTime.c
    Interpolation.F90
    MainUtils.F90
    Adaptive.F90
    EigenSolve.F90
    HashTable.F90
    MeshUtils.F90
    SParIterGlobals.F90
    SParIterComm.F90
    SParIterPrecond.F90
    SParIterSolver.F90
    Messages.F90
    Multigrid.F90
    Smoothers.F90
    ClusteringMethods.F90
    ParallelUtils.F90
    ParallelEigenSolve.F90
    solve_cmplx.F90
    solve_real.F90
    MGPrec.F90
    DefUtils.F90
    MeshPartition.F90
    MeshRemeshing.F90
    LUDecomposition.F90
    RadiationFactors.F90
    f_stubs.c
    ExchangeCorrelations.F90
    SolveHypre.c
    SolverActivate_x.F90
    SolveTrilinos.cxx
    SolveSuperLU.c
    iso_varying_string.F90
    umf4_f77wrapper.c
    VankaCreate.F90
    ParticleUtils.F90
    Feti.F90
    cholmod.c
    InterpolateMeshToMesh.F90
    InterpVarToVar.F90
    LinearForms.F90
    H1Basis.F90
    CircuitUtils.F90
    BackwardError.F90
    ElmerSolver.F90
    MagnetoDynamicsUtils.F90
    ComponentUtils.F90
    ZirkaHysteresis.F90)

set_property(SOURCE MaxwellAxiS.F90 PROPERTY COMPILE_DEFINITIONS FULL_INDUCTION)

if(WITH_LUA)
  include_directories(${LUA_INCLUDE_DIR})
  list(APPEND solverlib_SOURCES Lua.F90 elmer_lua_iface.c)
  install(
    DIRECTORY lua-scripts
    DESTINATION share/elmersolver/
    PATTERN lua-scripts/*.lua)
  file(GLOB DEFAULT_LUA_SCRIPTS "lua-scripts/*.lua")
  file(COPY ${DEFAULT_LUA_SCRIPTS}
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/share/elmersolver/lua-scripts/)
endif()

if(HAVE_FETI4I)
  list(APPEND solverlib_SOURCES ${FETI4I_INTERFACE_SOURCE})
endif()

# FILE(GLOB SRC_FILES *.src) FOREACH(FNAME ${SRC_FILES})
# GET_FILENAME_COMPONENT(BASENAME ${FNAME} NAME_WE) ADD_CUSTOM_COMMAND( OUTPUT
# ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.F90 COMMAND ${CMAKE_COMMAND} -E copy
# ${FNAME} ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.F90 DEPENDS
# ${CMAKE_CURRENT_SOURCE_DIR}/${BASENAME}.src ) ENDFOREACH()

include_directories("${CMAKE_BINARY_DIR}/fem/src")
include_directories("${CMAKE_SOURCE_DIR}/fhutiter/src")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/binio")
include_directories("${PROJECT_BINARY_DIR}/fhutiter/src")

if(WITH_Trilinos)
  include_directories("${Epetra_INCLUDE_DIRS}")
  include_directories("${ML_INCLUDE_DIRS}")
endif()

# Extract additional compile flags
get_directory_property(ELMER_CDEFS COMPILE_DEFINITIONS)
foreach(d ${ELMER_CDEFS})
  set(ELMER_F90FLAGS "${ELMER_F90FLAGS} -D${d}")
endforeach()

# Extract compile flags arising from build type
if(CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPCASE)
  set(ELMER_F90FLAGS
      "${ELMER_F90FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
  mark_as_advanced(CMAKE_BUILD_TYPE_UPCASE)
endif(CMAKE_BUILD_TYPE)

# Output elmerf90 and elmerld
set(HAVE_ELMERICE "FALSE")
if(WITH_ElmerIce)
  set(HAVE_ELMERICE "TRUE")
endif()

if(WIN32)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/elmerf90.bat.in.cmake
                 ${CMAKE_CURRENT_BINARY_DIR}/elmerf90.bat)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/elmerld.bat.in.cmake
                 ${CMAKE_CURRENT_BINARY_DIR}/elmerld.bat)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/elmerf90-nosh.in.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/elmerf90-nosh)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/elmerf90.in.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/elmerf90)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/elmerld.in.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/elmerld)

# Create pkg-config file if requested and place it under ${PKGCONFIG_PC_PATH}
# Defaults: PKGCONFIG_PC_PATH = ${CMAKE_INSTALL_PREFIX}/share/pkgconfig
if(CREATE_PKGCONFIG_FILE)
  set(PKGCONFIG_PC_PATH
      "${CMAKE_INSTALL_PREFIX}/share/pkgconfig"
      CACHE PATH "Path where to install elmer pkg-config file")
  configure_file(elmer.pc.cmake elmer.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/elmer.pc
          DESTINATION ${PKGCONFIG_PC_PATH})
endif(CREATE_PKGCONFIG_FILE)

# Copy elements.def and SOLVER.KEYWORDS (to enable testing)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/elements.def
               ${CMAKE_CURRENT_BINARY_DIR}/elements.def COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SOLVER.KEYWORDS
               ${CMAKE_CURRENT_BINARY_DIR}/SOLVER.KEYWORDS COPYONLY)

if(NOT MPI_FOUND)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mpif_stub.h
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/mpif_stub.h
       ${CMAKE_CURRENT_BINARY_DIR}/mpif.h)
endif()

add_library(mpi_stubs SHARED mpif_stubs.F90)
add_library(elmersolver SHARED ${solverlib_SOURCES})
set_target_properties(mpi_stubs PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(elmersolver PROPERTIES LINKER_LANGUAGE Fortran)

# We want MMG & Zoltan in the RPATH (if found), but given how many other
# components use ELMERSOLVER_RPATH_STRING, it seemed best to make a local copy
# just for ElmerSolver and libelmersolver
set(ELMERSOLVER_RPATH_STRING_MOD "${ELMERSOLVER_RPATH_STRING}")
if(HAVE_MMG)
  target_link_libraries(elmersolver ${MMG_LIBRARY})
  set(ELMERSOLVER_RPATH_STRING_MOD
      "${ELMERSOLVER_RPATH_STRING_MOD}:${MMG_LIBDIR}")
  set(ELMERLIB_RPATH_STRING "${ELMERLIB_RPATH_STRING}:${MMG_LIBDIR}")
  # TODO - add to RPath
endif()

if(HAVE_ZOLTAN)
  target_link_libraries(elmersolver ${ZOLTAN_LIBRARY})
  add_dependencies(elmersolver ${ZOLTAN_LIBRARY})
  # SET(ELMERSOLVER_RPATH_STRING_MOD
  # "${ELMERSOLVER_RPATH_STRING_MOD}/:${ZOLTAN_LIBDIR}")
  # SET(ELMERLIB_RPATH_STRING "${ELMERLIB_RPATH_STRING}/:${ZOLTAN_LIBDIR}")
endif()

# ElmerSolver libraries
set(ELMERSOLVER_LIBRARIES
    matc
    umfpack
    amd
    fhuti
    binio
    arpack
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${CMAKE_DL_LIBS})

if(WITH_LUA)
  list(APPEND ELMERSOLVER_LIBRARIES ${LUA_LIBRARIES})
endif()

add_executable(Solver_TGT Solver.F90)
set_target_properties(Solver_TGT PROPERTIES OUTPUT_NAME
                                            "${ELMERSOLVER_OUTPUT_NAME}")

if(NOT (WIN32))
  message(STATUS " ELMERSOLVER_RPATH_STRING_MOD "
                 ${ELMERSOLVER_RPATH_STRING_MOD})
  message(STATUS " ELMERLIB_RPATH_STRING " ${ELMERLIB_RPATH_STRING})

  set_target_properties(Solver_TGT PROPERTIES INSTALL_RPATH
                                              "${ELMERSOLVER_RPATH_STRING_MOD}")
  set_target_properties(elmersolver PROPERTIES INSTALL_RPATH
                                               "${ELMERLIB_RPATH_STRING}")
endif()

if(NOT (WITH_MPI))
  target_link_libraries(elmersolver mpi_stubs ${ELMERSOLVER_LIBRARIES})
  target_link_libraries(Solver_TGT elmersolver)
endif()

# ElmerSolver libraries
if(WITH_MPI)
  # ADD_LIBRARY(elmersolver_mpi SHARED ${solverlib_SOURCES})

  # Add parpack and possibly others
  if(Mumps_LIBRARIES)
    include_directories(${Mumps_INCLUDE_DIR})
    # Link to SCALAPACK implicitly
    list(APPEND ELMERSOLVER_LIBRARIES ${Mumps_LIBRARIES})
  elseif(MKL_CPARDISO_FOUND)
    # Mumps includes SCALAPACK, but it is needed with CPARDISO
    list(APPEND ELMERSOLVER_LIBRARIES ${SCALAPACK_LIBRARIES})
  endif()

  if(Hypre_LIBRARIES)
    list(APPEND ELMERSOLVER_LIBRARIES ${Hypre_LIBRARIES})
  endif()

  if(FETI4I_LIBRARIES)
    list(APPEND ELMERSOLVER_LIBRARIES ${FETI4I_LIBRARIES})
  endif()

  # Executable
  get_target_property(CURRENT_LINK_FLAGS Solver_TGT LINK_FLAGS)
  get_target_property(CURRENT_COMPILE_FLAGS Solver_TGT COMPILE_FLAGS)

  if(CURRENT_LINK_FLAGS)
    set(CURR_LFLAGS "${CURRENT_LINK_FLAGS};${MPI_Fortran_LINK_FLAGS}")
  else()
    set(CURR_LFLAGS "${MPI_Fortran_LINK_FLAGS}")
  endif()

  if(CURRENT_COMPILE_FLAGS)
    set(CURR_CFLAGS "${CURRENT_COMPILE_FLAGS};${MPI_Fortran_COMPILE_FLAGS}")
  else()
    set(CURR_CFLAGS "${MPI_Fortran_COMPILE_FLAGS}")
  endif()

  set_target_properties(Solver_TGT PROPERTIES LINK_FLAGS "${CURR_LFLAGS}")
  set_target_properties(Solver_TGT PROPERTIES COMPILE_FLAGS "${CURR_CFLAGS}")
  target_link_libraries(Solver_TGT elmersolver ${ELMERSOLVER_LIBRARIES})

  # Library object
  get_target_property(CURRENT_LINK_FLAGS elmersolver LINK_FLAGS)
  get_target_property(CURRENT_COMPILE_FLAGS elmersolver COMPILE_FLAGS)

  if(CURRENT_LINK_FLAGS)
    set(CURR_LFLAGS "${CURRENT_LINK_FLAGS};${MPI_Fortran_LINK_FLAGS}")
  else()
    set(CURR_LFLAGS "${MPI_Fortran_LINK_FLAGS}")
  endif()
  if(CURRENT_COMPILE_FLAGS)
    set(CURR_CFLAGS "${CURRENT_COMPILE_FLAGS};${MPI_Fortran_COMPILE_FLAGS}")
  else()
    set(CURR_CFLAGS "${MPI_Fortran_COMPILE_FLAGS}")
  endif()

  set_target_properties(elmersolver PROPERTIES LINK_FLAGS "${CURR_LFLAGS}")
  set_target_properties(elmersolver PROPERTIES COMPILE_FLAGS "${CURR_CFLAGS}")
  target_link_libraries(elmersolver ${ELMERSOLVER_LIBRARIES} parpack
                        ${MPI_Fortran_LIBRARIES})

  set(SERIAL_SOLVER_FILEPATH "${CMAKE_INSTALL_PREFIX}/bin/ElmerSolver")
  set(MPI_SOLVER_FILEPATH
      "${CMAKE_INSTALL_PREFIX}/bin/${ELMERSOLVER_OUTPUT_NAME}")
endif()

if(WITH_Trilinos)
  target_link_libraries(
    elmersolver
    "${Trilinos_LIBRARIES};${ML_LIBRARIES};${Epetra_LIBRARIES};${Teuchos_LIBRARIES}"
  )
endif()

# elmersolver -library TARGET_LINK_LIBRARIES(elmersolver
# ${ELMERSOLVER_LIBRARIES})

# ElmerSolver serial IF(NOT(NO_SERIAL_BINARY)) TARGET_LINK_LIBRARIES(elmersolver
# mpi_stubs ${ELMERSOLVER_LIBRARIES}) TARGET_LINK_LIBRARIES(ElmerSolver
# elmersolver) ENDIF()

# ElmerSolver ADD_EXECUTABLE(ElmerSolver Solver.F90)
# TARGET_LINK_LIBRARIES(ElmerSolver ${ELMERSOLVER_LIBRARIES} mpi_stubs)

# ViewFactors
add_executable(ViewFactors ViewFactors.F90)
set_target_properties(ViewFactors PROPERTIES LINKER_LANGUAGE Fortran)
# TARGET_INCLUDE_DIRECTORIES(ViewFactors PUBLIC
# ${CMAKE_BINARY_DIR}/fem/src/view3d) TARGET_INCLUDE_DIRECTORIES(ViewFactors
# PUBLIC ${CMAKE_BINARY_DIR}/fem/src/viewfactors) FIND_LIBRARY(VIEW3D_LIBRARIES
# PATHS ${CMAKE_CURRENT_BINARY_DIR}/view3d) FIND_LIBRARY(VIEWAXIS_LIBRARIES
# PATHS ${CMAKE_CURRENT_BINARY_DIR}/viewaxis)
target_link_libraries(ViewFactors ${ELMERSOLVER_LIBRARIES} mpi_stubs view3d
                      viewaxis elmersolver)
install(TARGETS ViewFactors RUNTIME DESTINATION "bin")

if(NOT (WIN32))
  set_target_properties(ViewFactors PROPERTIES INSTALL_RPATH
                                               "${ELMERSOLVER_RPATH_STRING}")
endif()

# GebhardFactors
add_executable(GebhardtFactors GebhardtFactors.F90)
target_link_libraries(GebhardtFactors ${ELMERSOLVER_LIBRARIES} mpi_stubs
                      elmersolver)

if(NOT (WIN32))
  set_target_properties(GebhardtFactors
                        PROPERTIES INSTALL_RPATH "${ELMERSOLVER_RPATH_STRING}")
endif()

if(WITH_Trilinos)
  target_link_libraries(
    GebhardtFactors
    elmersolver
    "${Trilinos_LIBRARIES};${Belos_LIBRARIES};${ML_LIBRARIES};${Epetra_LIBRARIES};${Teuchos_LIBRARIES}"
  )
endif()

install(
  TARGETS GebhardtFactors mpi_stubs
  RUNTIME DESTINATION "bin"
  LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}")

# Installation rules
install(FILES elements.def SOLVER.KEYWORDS DESTINATION "share/elmersolver/lib")

if(NOT (WIN32))
  if(WITH_MPI)
    install(
      TARGETS Solver_TGT elmersolver
      RUNTIME DESTINATION "bin"
      LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}")
    install(
      TARGETS Solver_TGT
      RUNTIME DESTINATION "bin"
      LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}")
  else()
    install(
      TARGETS Solver_TGT elmersolver mpi_stubs
      RUNTIME DESTINATION "bin"
      LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}")
  endif()
else()
  install(
    TARGETS elmersolver mpi_stubs
    RUNTIME DESTINATION "${ELMER_INSTALL_LIB_DIR}"
    LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}")
  install(
    TARGETS Solver_TGT elmersolver mpi_stubs
    RUNTIME DESTINATION "bin"
    LIBRARY DESTINATION "${ELMER_INSTALL_LIB_DIR}")
endif()

if(WIN32)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/elmerf90.bat DESTINATION "bin")
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/elmerld.bat DESTINATION "bin")
endif()
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/elmerf90 DESTINATION "bin")
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/elmerld DESTINATION "bin")

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/.
        DESTINATION "share/elmersolver/include")

if(WITH_MPI)
  if(NOT (WIN32))
    # INSTALL(CODE " EXECUTE_PROCESS(COMMAND \${CMAKE_COMMAND} -E create_symlink
    # ElmerSolver_mpi ElmerSolver WORKING_DIRECTORY
    # \${CMAKE_INSTALL_PREFIX}/bin)")
    add_custom_target(
      elmersolver_link ALL COMMAND ${CMAKE_COMMAND} -E create_symlink
                                   ElmerSolver_mpi ElmerSolver)
    add_dependencies(elmersolver_link Solver_TGT)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ElmerSolver DESTINATION "bin")
  else()
    install(CODE "EXECUTE_PROCESS(COMMAND \${CMAKE_COMMAND} -E copy
    \${CMAKE_INSTALL_PREFIX}/bin/ElmerSolver_mpi.exe
    \${CMAKE_INSTALL_PREFIX}/bin/ElmerSolver.exe)")
  endif()
endif()
