INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( cyl-case.sif cyl-case.sif COPYONLY)
ADD_ELMERTEST_MODULE(mgdyn_anisotropic_rel reluctivity reluctivity.F90)

file(COPY ELMERSOLVER_STARTINFO reluctivity.F90 cylinders DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_anisotropic_rel LABELS whitney)
