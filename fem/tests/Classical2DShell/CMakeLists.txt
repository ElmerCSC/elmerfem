INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( cylindrical_shell.sif cylindrical_shell.sif COPYONLY)
ADD_ELMERTEST_MODULE(Classical2DShell BenchmarkShell BenchmarkShell.F90)

file(COPY ELMERSOLVER_STARTINFO eighth.grd BenchmarkShell.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Classical2DShell LABELS quick)
