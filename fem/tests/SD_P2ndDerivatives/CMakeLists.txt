INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( ptest.sif ptest.sif COPYONLY)
ADD_ELMERTEST_MODULE(SD_P2ndDerivatives P2ndDerivativesTest P2ndDerivativesTest.f90)

file(COPY P2ndDerivativesTest.f90 ELMERSOLVER_STARTINFO mesh  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(SD_P2ndDerivatives)
ADD_ELMER_LABEL(SD_P2ndDerivatives quick serendipity)
