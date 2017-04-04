INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( Bentonite.sif Bentonite.sif COPYONLY)
ADD_ELMERTEST_MODULE(bentonite BentoniteSolver BentoniteSolver.F90)

file(COPY ELMERSOLVER_STARTINFO  BentoniteSolver.F90 1d DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(bentonite)
ADD_ELMER_LABEL(bentonite quick)
