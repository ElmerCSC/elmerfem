INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(SD_BlockLinElast1 BlockLinElast BlockLinElast.f90)

file(COPY square.grd ELMERSOLVER_STARTINFO BlockLinElast.f90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(SD_BlockLinElast1 LABELS quick block serendipity)
