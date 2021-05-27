include(test_macros)
file(MAKE_DIRECTORY Ldomain)
execute_process(COMMAND ${MESH2D_BIN} Ldomain.mif .)
file(COPY mesh.boundary mesh.elements mesh.header mesh.nodes
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Ldomain)
file(REMOVE mesh.boundary mesh.elements mesh.header mesh.nodes)
file(COPY Ldomain.mif DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Ldomain)
RUN_ELMER_TEST()
