INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(KeywordCompareTensor KeywordCompareTensor KeywordCompareTensor.F90)

file(COPY ELMERSOLVER_STARTINFO KeywordCompareTensor.F90 square.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(KeywordCompare5 LABELS matc)
