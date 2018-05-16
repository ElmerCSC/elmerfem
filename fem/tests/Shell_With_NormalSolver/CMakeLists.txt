INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( shell_case1.sif shell_case1.sif COPYONLY)
ADD_ELMERTEST_MODULE(Shell_With_NormalSolver CylinderNormalSolver CylinderNormalSolver.F90)

file(COPY ELMERSOLVER_STARTINFO p1mesh3 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Shell_With_NormalSolver NPROCS 1 2 LABELS shell)
