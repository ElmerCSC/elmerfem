INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(mgdyn_harmonic_wire_Cgauge_automatic)


file(COPY ELMERSOLVER_STARTINFO wire.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_harmonic_wire_Cgauge_automatic LABELS whitney vtu)
