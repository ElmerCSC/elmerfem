INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)
CONFIGURE_FILE(sif/2241.sif sif/2241.sif COPYONLY) 
file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY sif/2241-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
file(COPY 2241/mesh.boundary 2241/mesh.header 2241/mesh.nodes 2241/mesh.elements 2241/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241")
ADD_ELMER_TEST(circuits_harmonic_stranded_homogenization LABELS 3D circuits harmonic homogenization stranded mgdyn whitney)

