INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( waveguide.sif waveguide.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO shoebox_tets DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(VectorHelmholtzWaveguide)
ADD_ELMER_LABEL(VectorHelmholtzWaveguide benchmark)
