include(${TEST_SOURCE}/../test_macros.cmake)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 mortar.msh -autoclean)
execute_process(COMMAND ${ELMERGRID_BIN} extrude.eg)

RUN_ELMER_TEST_NEW()
