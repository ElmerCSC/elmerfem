INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( elasticity.sif elasticity.sif COPYONLY)

file(COPY elasticity.grd ELMERSOLVER_STARTINFO  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(elasticity)
ADD_ELMER_LABEL(elasticity quick)
