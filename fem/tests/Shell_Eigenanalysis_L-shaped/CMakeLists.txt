INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( shell_eigenanalysis.sif shell_eigenanalysis.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO L DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Shell_Eigenanalysis_L-shaped LABELS shell eigen)
