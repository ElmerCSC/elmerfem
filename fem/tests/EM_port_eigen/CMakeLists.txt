INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( port_eigenanalysis.sif port_eigenanalysis.sif COPYONLY)

file(COPY port.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(EM_port_eigen LABELS eigen complex_eigen)
