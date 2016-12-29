INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( zalesak.sif zalesak.sif COPYONLY)
ADD_ELMERTEST_MODULE(ParticleAdvectorZalesak Zalesak Zalesak.f90)

file(COPY square.grd ELMERSOLVER_STARTINFO Zalesak.f90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ParticleAdvectorZalesak)
ADD_ELMER_LABEL(ParticleAdvectorZalesak particle)
