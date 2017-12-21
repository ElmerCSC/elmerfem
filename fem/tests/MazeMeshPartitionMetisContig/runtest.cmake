
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 maze.grd -metis 6 -metiscontig -nooverwrite)

RUN_ELMER_TEST()
