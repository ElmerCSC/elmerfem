INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(RigidMeshMapperCatenary CatenaryCurve CatenaryCurve.F90)

file(COPY ELMERSOLVER_STARTINFO CatenaryCurve.F90 rectangles.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(RigidMeshMapperCatenary)
