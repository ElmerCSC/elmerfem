include(${TEST_SOURCE}/../test_macros.cmake)

execute_process(COMMAND 
  ${CMAKE_COMMAND} -E copy_directory 
  ${TEST_SOURCE}/1d ${CMAKE_CURRENT_BINARY_DIR}/1d)

RUN_ELMER_TEST()
