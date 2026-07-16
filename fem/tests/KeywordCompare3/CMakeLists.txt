INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(KeywordCompareRealVec KeywordCompareRealVec KeywordCompareRealVec.F90)
ADD_ELMERTEST_MODULE(RealVecFunctionTests RealVecFunctionTests RealVecFunctionTests.F90)

file(COPY ELMERSOLVER_STARTINFO KeywordCompareRealVec.F90 RealVecFunctionTests.F90 square.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(KeywordCompare3 LABELS matc)
