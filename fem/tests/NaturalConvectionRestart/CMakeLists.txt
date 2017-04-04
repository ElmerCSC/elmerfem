INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
CONFIGURE_FILE( init.sif init.sif COPYONLY)

file(COPY square.grd init_asc.result ELMERSOLVER_STARTINFO ! DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(NaturalConvectionRestart)
