INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(h1basistest.sif h1basistest.sif COPYONLY)
ADD_ELMERTEST_MODULE(H1BasisEvaluation H1BasisEvaluation H1BasisEvaluation.F90)

file(COPY H1BasisEvaluation.F90 cube.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(H1BasisEvaluation LABELS benchmark serial)
