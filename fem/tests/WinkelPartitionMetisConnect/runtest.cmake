include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 winkel.grd -partdual -metisrec 8 -connect 3 7 -partlayers 3 -nooverwrite)

RUN_ELMER_TEST()
