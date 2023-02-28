INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( test.sif test.sif COPYONLY)
ADD_ELMERTEST_MODULE(L2ndDerivatives L2ndDerivativesTest L2ndDerivativesTest.f90)

file(COPY L2ndDerivativesTest.f90 ELMERSOLVER_STARTINFO mesh  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(L2ndDerivatives)
ADD_ELMER_LABEL(L2ndDerivatives quick)
