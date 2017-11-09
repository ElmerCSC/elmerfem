INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(ParticleFallingBlock CalculateFraction CalculateFraction.F90)

file(COPY block.grd CalculateFraction.F90 ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ParticleFallingBlock)
ADD_ELMER_LABEL(ParticleFallingBlock particle)
