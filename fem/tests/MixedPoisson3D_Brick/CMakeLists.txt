INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( mixed.sif mixed.sif COPYONLY)
ADD_ELMERTEST_MODULE(MixedPoisson3D_Brick Source Source.F90)

file(COPY ELMERSOLVER_STARTINFO testcube.grd Source.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(MixedPoisson3D_Brick LABELS quick vector_element)
