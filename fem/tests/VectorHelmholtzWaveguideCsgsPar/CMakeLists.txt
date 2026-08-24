INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(waveguide.sif waveguide.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO cjisgs.sif shoebox_prisms.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(VectorHelmholtzWaveguideCsgsPar NPROCS 1 2 3 LABELS whitney emwave)
