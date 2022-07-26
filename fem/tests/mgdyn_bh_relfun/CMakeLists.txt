INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(mgdyn_bh_relfun reluctivity reluctivity.F90)

file(COPY ELMERSOLVER_STARTINFO mesh_tet HB.dat reluctivity.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_bh_relfun LABELS whitney)
