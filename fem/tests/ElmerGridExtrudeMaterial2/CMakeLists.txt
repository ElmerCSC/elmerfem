INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO squares.grd squares2cubes.eg DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ElmerGridExtrudeMaterial2)
ADD_ELMER_LABEL(ElmerGridExtrudeMaterial2 quick)
