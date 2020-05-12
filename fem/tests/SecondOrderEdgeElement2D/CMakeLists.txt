INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( convergencetest.sif convergencetest.sif COPYONLY)
ADD_ELMERTEST_MODULE(SecondOrderEdgeElement2D EdgeElementTest EdgeElementTest.F90)

file(COPY ELMERSOLVER_STARTINFO square.grd EdgeElementTest.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(SecondOrderEdgeElement2D LABELS quick vector_element)
