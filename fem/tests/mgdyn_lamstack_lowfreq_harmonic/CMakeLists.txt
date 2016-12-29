INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( lamstack_lowfreq.sif lamstack_lowfreq.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO toroidalcoil DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_lamstack_lowfreq_harmonic LABELS whitney)
