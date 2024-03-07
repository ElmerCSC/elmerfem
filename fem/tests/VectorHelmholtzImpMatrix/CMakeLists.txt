INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case4.sif case4.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO four_channels.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(VectorHelmholtzImpMatrix LABELS whitney cmodes lumping)
