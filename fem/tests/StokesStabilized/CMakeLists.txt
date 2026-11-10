INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(quad.sif quad.sif COPYONLY)
CONFIGURE_FILE(tri.sif  tri.sif  COPYONLY)
CONFIGURE_FILE(tet.sif  tet.sif  COPYONLY)
CONFIGURE_FILE(quad_taylorhood.sif quad_taylorhood.sif COPYONLY)
CONFIGURE_FILE(tri_taylorhood.sif  tri_taylorhood.sif  COPYONLY)
CONFIGURE_FILE(tet_taylorhood.sif  tet_taylorhood.sif  COPYONLY)

file(COPY Step.grd tet.msh ELMERSOLVER_STARTINFO
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(StokesStabilized LABELS quick)
