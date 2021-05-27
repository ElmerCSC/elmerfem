include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 FlatPlate.msh -autoclean)
RUN_ELMER_TEST(ELMER_LIB contrib)
