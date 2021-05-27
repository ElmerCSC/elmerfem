INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( mixed_nedelec2nd.sif mixed_nedelec2nd.sif COPYONLY)
ADD_ELMERTEST_MODULE(MixedPoisson3D Source Source.F90)

file(COPY ELMERSOLVER_STARTINFO cube Source.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(MixedPoisson3D LABELS quick vector_element)
