IF(HAVE_MMG AND HAVE_ZOLTAN)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

  CONFIGURE_FILE(calvingMMG.sif calvingMMG.sif COPYONLY)
  CONFIGURE_FILE(LeftRail.xy LeftRail.xy COPYONLY)
  CONFIGURE_FILE(RightRail.xy RightRail.xy COPYONLY)

  ADD_ELMERICETEST_MODULE(Calving3D_lset bedrockfunction_3d ./PROG/bedrockfunction_3D.F90)
  ADD_ELMERICETEST_MODULE(Calving3D_lset DummySolver ./PROG/DummySolver.F90)

  FILE(COPY ELMERSOLVER_STARTINFO PlanMesh.msh PROG DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Results/PlanMesh")
  FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Results/isoline")
  FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Results/calving_plane")

  IF(WITH_MPI)
    # Check the task bounds and add only compatible tests
    SET(NPROCS 2)
    IF(${NPROCS} GREATER ${MPI_TEST_MAXPROC} OR ${NPROCS} LESS ${MPI_TEST_MINPROC})
      MESSAGE(STATUS "Skipping test Calving3D_lset with ${NPROCS} procs")
    ELSE()
      ADD_ELMERICE_TEST(Calving3D_lset LABELS slow)
    ENDIF()
  ENDIF()

ENDIF()
