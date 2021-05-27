INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(SphereDistance SphereDistance SphereDistance.F90)

file(COPY ELMERSOLVER_STARTINFO SphereDistance.F90 cube.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(MeshSolverCubeToSphere)
