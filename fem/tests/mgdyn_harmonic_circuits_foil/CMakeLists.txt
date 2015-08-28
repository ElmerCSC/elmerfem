INCLUDE(${CMAKE_SOURCE_DIR}/fem/tests/test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)
CONFIGURE_FILE(sif/432.sif sif/432.sif COPYONLY) 
file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY sif/432-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
file(COPY 432/mesh.boundary 432/mesh.header 432/mesh.nodes 432/mesh.elements 432/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/432")
ADD_ELMER_TEST(circuits3D_harmonic_foil)
ADD_ELMER_LABEL(circuits3D_harmonic_foil circuits3D_harmonic)

