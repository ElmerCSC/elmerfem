include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step)
file(COPY flow.dat flow.dat.pos DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Step/)
RUN_ELMER_TEST()
