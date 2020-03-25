INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( with_springs.sif with_springs.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO beam DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(elasticity_with_springs LABELS quick elasticsolve)
