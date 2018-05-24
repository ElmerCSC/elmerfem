IF(WITH_Hypre)
  INCLUDE(test_macros)

  CONFIGURE_FILE(case.sif case.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO cylinders DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_hypre_ams NPROCS 2 LABELS quick hypre)
ENDIF()
