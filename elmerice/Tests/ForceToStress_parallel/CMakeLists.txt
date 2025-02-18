INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(test.sif test.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO cube.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

IF(WITH_MPI)
  # Check the task bounds and add only compatible tests
  SET(NPROCS 4)
  IF(${NPROCS} GREATER ${MPI_TEST_MAXPROC} OR ${NPROCS} LESS ${MPI_TEST_MINPROC})
    MESSAGE(STATUS "Skipping test ForceToStress_parallel with ${NPROCS} procs")
  ELSE()
    ADD_ELMERICE_TEST(ForceToStress_parallel)
    ADD_ELMERICE_LABEL(ForceToStress_parallel elmerice-fast)
  ENDIF()
ENDIF()
