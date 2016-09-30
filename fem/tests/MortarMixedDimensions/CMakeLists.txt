INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)

file(COPY square.grd ELMERSOLVER_STARTINFO cube.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(MortarMixedDimensions)
ADD_ELMER_LABEL(MortarMixedDimensions quick)
