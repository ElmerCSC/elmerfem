INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( constant_pressure_load.sif constant_pressure_load.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO Mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Beam_3D_Cantilever_Orientation LABELS beam)
