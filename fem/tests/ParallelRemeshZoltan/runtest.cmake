INCLUDE(test_macros)

SET(NPROCS 4)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 cube.msh -autoclean -partdual -metiskway ${NPROCS})

RUN_ELMER_TEST()
