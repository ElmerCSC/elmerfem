INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO extrude_quarter.eg polygons.m mortar_quarter.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(RotatingBCPoisson3DGeneric)
ADD_ELMER_LABEL(RotatingBCPoisson3DGeneric benchmark)
