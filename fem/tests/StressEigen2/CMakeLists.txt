INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO angle DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(StressEigen2 NPROCS 1 4 8 LABELS quick eigen)

# Tag the cases which run subroutines from ParallelEigenSolve.F90
FOREACH(_np 4 8)
  IF(TEST StressEigen2_np${_np})
    SET_PROPERTY(TEST StressEigen2_np${_np} APPEND PROPERTY LABELS eigen_parallel)
  ENDIF()
ENDFOREACH()
