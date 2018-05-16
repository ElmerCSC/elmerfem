INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(KeywordCompareReal KeywordCompareReal KeywordCompareReal.F90)
ADD_ELMERTEST_MODULE(RealFunctionTests RealFunctionTests RealFunctionTests.F90)

file(COPY ELMERSOLVER_STARTINFO KeywordCompareReal.F90 RealFunctionTests.F90 square.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(KeywordCompare2 LABELS matc)
