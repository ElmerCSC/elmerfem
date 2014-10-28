INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( levy.sif levy.sif COPYONLY)

file(COPY mesh.elements mesh.nodes ELMERSOLVER_STARTINFO mesh.header levy.grd mesh.boundary DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(StrainCalculation02)
