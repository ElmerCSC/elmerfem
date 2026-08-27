include(test_macros)

execute_process(COMMAND ${ELMERGRID_BIN} 1 2 coil -relh 1.0 -out coarse)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 coil -relh 0.6 -out fine)

RUN_ELMER_TEST()
