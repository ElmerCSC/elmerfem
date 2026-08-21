INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

# The checker, the profiles and the mesh are taken from the sibling test rather
# than copied, so the control cannot drift away from the case it controls for.
SET(MMS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../StatCurrentDiffusionCurrent)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(StatCurrentDiffusionCurrentOff Profiles ${MMS_SOURCE_DIR}/Profiles.F90)
ADD_ELMERTEST_MODULE(StatCurrentDiffusionCurrentOff MMSCheck ${MMS_SOURCE_DIR}/MMSCheck.F90)

file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY ${MMS_SOURCE_DIR}/strip.grd
          ${MMS_SOURCE_DIR}/Profiles.F90
          ${MMS_SOURCE_DIR}/MMSCheck.F90
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(StatCurrentDiffusionCurrentOff LABELS quick electrochemistry)
