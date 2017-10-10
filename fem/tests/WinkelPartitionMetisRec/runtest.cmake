execute_process(COMMAND ${ELMERGRID_BIN} 1 2 winkel.grd -partdual -metisrec 8 -nooverwrite)

RUN_ELMER_TEST()
