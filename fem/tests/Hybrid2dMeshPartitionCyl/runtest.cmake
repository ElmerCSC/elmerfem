include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 circle_in_box.grd -partcyl 1 8 1 -nooverwrite)
RUN_ELMER_TEST()
