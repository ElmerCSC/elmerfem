INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( cavity_parstokes.sif cavity_parstokes.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO sinus0.grd  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ParStokes_NormalTangentialBC)
ADD_ELMER_LABEL(ParStokes_NormalTangentialBC quick)
