INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( 3dcavity_parstokes.sif 3dcavity_parstokes.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO sinus0.grd extrude.eg DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ParStokes_NormalTangentialBC3D)

