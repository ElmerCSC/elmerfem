INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(fsi_beam_optimize WallProfiles WallProfiles.f90)

file(COPY beam_channel.in2d ELMERSOLVER_STARTINFO case_b.best0 WallProfiles.f90 mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(fsi_beam_optimize LABELS elasticsolve)
