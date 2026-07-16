IF(HAVE_XIOS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

  CONFIGURE_FILE(Case.sif Case.sif COPYONLY)

  ADD_ELMERICETEST_MODULE(Xios2 Check_xios2 Check_xios2.F90)
  FILE(COPY Check_xios2.F90 rectangle.grd ELMERSOLVER_STARTINFO iodef.xml context_elmerice.xml  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  IF(WITH_MPI)
    # Check the task bounds and add only compatible tests
    SET(NPROCS 2)
    IF(${NPROCS} GREATER ${MPI_TEST_MAXPROC} OR ${NPROCS} LESS ${MPI_TEST_MINPROC})
      MESSAGE(STATUS "Skipping test Xios2 with ${NPROCS} procs")
    ELSE()
      ADD_ELMERICE_TEST(Xios2)
      ADD_ELMERICE_LABEL(Xios2 elmerice-fast)
      ADD_ELMERICE_LABEL(Xios2 elmerice-xios)
      ADD_ELMERICE_LABEL(Xios2 netcdf)
    ENDIF()
  ENDIF()
ENDIF()
