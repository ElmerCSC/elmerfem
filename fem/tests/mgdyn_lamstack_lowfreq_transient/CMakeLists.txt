INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( lamstack.sif lamstack.sif COPYONLY)
ADD_ELMERTEST_MODULE(mgdyn_lamstack_lowfreq_transient currents currents.f90)

file(COPY ELMERSOLVER_STARTINFO currents.f90 toroidalcoil DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_lamstack_lowfreq_transient LABELS whitney)
