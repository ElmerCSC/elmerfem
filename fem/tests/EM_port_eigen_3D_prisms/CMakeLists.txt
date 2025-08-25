INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( port_3D_nodalprec.sif port_3D_nodalprec.sif COPYONLY)

file(COPY port.grd extrude.eg ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(EM_port_eigen_3D_prisms LABELS eigen complex_eigen)
