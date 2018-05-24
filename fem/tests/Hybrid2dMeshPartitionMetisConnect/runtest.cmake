execute_process(COMMAND ${ELMERGRID_BIN} 1 2 circle_in_box.grd -metisrec 8 -connect 2 -partlayers 3 -partdual -nooverwrite)

RUN_ELMER_TEST()
