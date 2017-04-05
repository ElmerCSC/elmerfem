INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( levy.sif levy.sif COPYONLY)

file(COPY mesh.elements mesh.nodes ELMERSOLVER_STARTINFO mesh.header levy.dat mesh.boundary DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(StrainCalculation01)
ADD_ELMER_LABEL(StrainCalculation01 quick)
