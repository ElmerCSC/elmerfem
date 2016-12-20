execute_process(COMMAND ${ELMERGRID_BIN} 14 2 sector.msh -autoclean)
execute_process(COMMAND ${ELMERGRID_BIN} sector_extrude.eg)

RUN_ELMER_TEST()
