INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( interpolationtest.sif interpolationtest.sif COPYONLY)
ADD_ELMERTEST_MODULE(SecondOrderEdgeElement2D_BCs EdgeFEInterpolationTest EdgeFEInterpolationTest.F90)

file(COPY ELMERSOLVER_STARTINFO square.grd EdgeFEInterpolationTest.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(SecondOrderEdgeElement2D_BCs)
ADD_ELMER_LABEL(SecondOrderEdgeElement2D_BCs quick)
