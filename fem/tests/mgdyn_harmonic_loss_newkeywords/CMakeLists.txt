INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( MGDynamicsHarmonic.sif MGDynamicsHarmonic.sif COPYONLY)

# The mesh is shared with the mgdyn_harmonic_loss test
file(COPY ELMERSOLVER_STARTINFO ${CMAKE_CURRENT_SOURCE_DIR}/../mgdyn_harmonic_loss/mesh3d DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_harmonic_loss_newkeywords LABELS whitney)
