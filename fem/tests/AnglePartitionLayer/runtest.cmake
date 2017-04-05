execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle.grd -partition 3 1 1 -connect 1 2 -nooverwrite)

RUN_ELMER_TEST()
