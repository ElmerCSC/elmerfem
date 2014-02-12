include(${TEST_SOURCE}/../test_macros.cmake)

execute_process(COMMAND 
  ${CMAKE_COMMAND} -E copy_directory 
  ${TEST_SOURCE}/Ldomain ${CMAKE_CURRENT_BINARY_DIR}/Ldomain)

RUN_ELMER_TEST()
