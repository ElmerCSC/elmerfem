INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(Stokes_prognostic_nonlWeert.sif Stokes_prognostic_nonlWeert.sif COPYONLY)
FILE(COPY ELMERSOLVER_STARTINFO README mesh2D perlin.lua bedrock.dat surf.dat DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(FSSA_perlin2d)
ADD_ELMERICE_LABEL(FSSA_perlin2d elmerice-long)
ADD_ELMERICE_LABEL(FSSA_perlin2d slow)
