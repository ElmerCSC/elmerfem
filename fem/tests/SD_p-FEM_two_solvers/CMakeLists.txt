INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( two_solvers.sif two_solvers.sif COPYONLY)
ADD_ELMERTEST_MODULE(SD_p-FEM_two_solvers BenchmarkShell BenchmarkShell.F90)

file(COPY ELMERSOLVER_STARTINFO eighth.grd BenchmarkShell.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(SD_p-FEM_two_solvers LABELS quick p-fem serendipity)
