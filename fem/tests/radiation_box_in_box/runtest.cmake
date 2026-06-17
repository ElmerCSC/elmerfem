include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 box_in_box.grd)
RUN_ELMER_TEST()
