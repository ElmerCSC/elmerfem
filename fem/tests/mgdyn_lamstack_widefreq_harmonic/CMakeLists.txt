INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( lamstack_widefreq.sif lamstack_widefreq.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO toroidalcoil DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_lamstack_widefreq_harmonic LABELS whitney)
