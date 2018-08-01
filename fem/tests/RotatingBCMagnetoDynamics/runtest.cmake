include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 mortar_quarter.msh -autoclean)
execute_process(COMMAND ${ELMERGRID_BIN} extrude_quarter.eg)

RUN_ELMER_TEST()
