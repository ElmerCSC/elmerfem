include(${TEST_SOURCE}/../test_macros.cmake)

execute_process(COMMAND
  ${MESH2D_BIN} Ldomain.mif Ldomain
  )
# execute_process(COMMAND 
#   ${CMAKE_COMMAND} -E copy 
#   ${CMAKE_CURRENT_BINARY_DIR}/mesh.* 
#   ${CMAKE_CURRENT_BINARY_DIR}/Ldomain)
execute_process(COMMAND 
  ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_BINARY_DIR}/Ldomain.mif 
  ${CMAKE_CURRENT_BINARY_DIR}/Ldomain)

RUN_ELMER_TEST()