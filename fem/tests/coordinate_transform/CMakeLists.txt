INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)
CONFIGURE_FILE(sif/1966.sif sif/1966.sif COPYONLY) 
file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY sif/1966-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
file(COPY 1966/mesh.boundary 1966/mesh.header 1966/mesh.nodes 1966/mesh.elements 1966/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1966")
ADD_ELMER_TEST(coordinate_transform)
ADD_ELMER_LABEL(coordinate_transform circuit_utils)

