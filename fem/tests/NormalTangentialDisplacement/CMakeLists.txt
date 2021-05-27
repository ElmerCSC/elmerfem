INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(NormalTangentialDisplacement TangentForce TangentForce.f90)

file(COPY square.grd ELMERSOLVER_STARTINFO TangentForce.f90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(NormalTangentialDisplacement LABELS quick elasticsolve)

