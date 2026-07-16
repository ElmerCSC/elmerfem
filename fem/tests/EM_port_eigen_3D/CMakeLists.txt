INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( port_eigenanalysis_3D.sif port_eigenanalysis_3D.sif COPYONLY)

file(COPY port.grd extrude.eg ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(EM_port_eigen_3D LABELS eigen complex_eigen)
