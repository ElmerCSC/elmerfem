IF(HAVE_XIOS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

  CONFIGURE_FILE(Case.sif Case.sif COPYONLY)

  ADD_ELMERICETEST_MODULE(Xios Check Check.F90)
  FILE(COPY Check.F90 rectangle.grd ELMERSOLVER_STARTINFO iodef.xml context_elmerice.xml  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  IF(WITH_MPI)
    # Check the task bounds and add only compatible tests
    SET(NPROCS 2)
    IF(${NPROCS} GREATER ${MPI_TEST_MAXPROC} OR ${NPROCS} LESS ${MPI_TEST_MINPROC})
      MESSAGE(STATUS "Skipping test Xios with ${NPROCS} procs")
    ELSE()
      ADD_ELMERICE_TEST(Xios)
      ADD_ELMERICE_LABEL(Xios elmerice-fast)
      ADD_ELMERICE_LABEL(Xios elmerice-xios)
      ADD_ELMERICE_LABEL(Xios netcdf)
    ENDIF()
  ENDIF()

ENDIF()
