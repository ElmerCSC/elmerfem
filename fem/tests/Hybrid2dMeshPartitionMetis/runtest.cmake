
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 circle_in_box.grd -metis 8 -partdual -nooverwrite)

RUN_ELMER_TEST()
