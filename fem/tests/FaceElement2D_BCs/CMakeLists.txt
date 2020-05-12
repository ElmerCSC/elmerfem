INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( interpolationtest.sif interpolationtest.sif COPYONLY)
ADD_ELMERTEST_MODULE(FaceElement2D_BCs FaceFEInterpolationTest FaceFEInterpolationTest.F90)

file(COPY ELMERSOLVER_STARTINFO square.grd FaceFEInterpolationTest.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(FaceElement2D_BCs LABELS quick vector_element)
