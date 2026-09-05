INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( waveguide.sif waveguide.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO shoebox_tets DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(VectorHelmholtzWaveguidePar NPROCS 2 4 LABELS whitney emwave complex)
