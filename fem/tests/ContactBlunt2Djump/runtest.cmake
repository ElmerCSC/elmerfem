execute_process(COMMAND ${ELMERGRID_BIN} 14 2 blunt.msh -autoclean )
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 base.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 blunt -in base -unite -out mesh)

RUN_ELMER_TEST()
