INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(StatCurrentDiffusionCurrent Profiles Profiles.F90)
ADD_ELMERTEST_MODULE(StatCurrentDiffusionCurrent MMSCheck MMSCheck.F90)

file(COPY ELMERSOLVER_STARTINFO strip.grd Profiles.F90 MMSCheck.F90
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(StatCurrentDiffusionCurrent LABELS quick electrochemistry)
