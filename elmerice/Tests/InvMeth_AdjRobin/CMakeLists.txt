INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(Robin_Beta.sif Robin_Beta.sif COPYONLY)

ADD_ELMERICETEST_MODULE(InvMeth_AdjRobin USF_Init ./PROG/USF_Init.f90)

FILE(COPY mesh2D.grd ELMERSOLVER_STARTINFO PROG DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

IF(WITH_MPI)
  # Check the task bounds and add only compatible tests
  SET(NPROCS 4)
  IF(${NPROCS} GREATER ${MPI_TEST_MAXPROC} OR ${NPROCS} LESS ${MPI_TEST_MINPROC})
    MESSAGE(STATUS "Skipping test InvMeth_AdjRobin with ${NPROCS} procs")
  ELSE()
    ADD_ELMERICE_TEST(InvMeth_AdjRobin)
    ADD_ELMERICE_LABEL(InvMeth_AdjRobin elmerice-long)
    ADD_ELMERICE_LABEL(InvMeth_AdjRobin slow)
  ENDIF()
ENDIF()
