execute_process(COMMAND ${ELMERGRID_BIN} 2 2 cylinders -metis 2 3 -nooverwrite)

RUN_ELMER_TEST()
