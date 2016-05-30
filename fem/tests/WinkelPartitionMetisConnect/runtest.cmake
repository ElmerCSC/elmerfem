execute_process(COMMAND ${ELMERGRID_BIN} 1 2 winkel.grd -partdual -metis 8 3 -connect 3 7 -nooverwrite)

RUN_ELMER_TEST()
