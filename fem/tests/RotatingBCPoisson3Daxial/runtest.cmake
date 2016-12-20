execute_process(COMMAND ${ELMERGRID_BIN} 14 2 sector.msh -autoclean)
execute_process(COMMAND ${ELMERGRID_BIN} extrude_sector.eg)

RUN_ELMER_TEST()
