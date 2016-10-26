INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( toroidalcoil.sif toroidalcoil.sif COPYONLY)
ADD_ELMERTEST_MODULE(mgdyn_anisotropic_cond currents currents.f90)

file(COPY ELMERSOLVER_STARTINFO currents.f90 toroidalcoil DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_anisotropic_cond LABELS whitney)
