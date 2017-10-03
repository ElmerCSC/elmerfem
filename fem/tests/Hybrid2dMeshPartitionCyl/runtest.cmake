execute_process(COMMAND ${ELMERGRID_BIN} 1 2 circle_in_box.grd -partition 1 8 1 3 -nooverwrite)
RUN_ELMER_TEST()
