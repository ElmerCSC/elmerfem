INCLUDE(test_macros)

SET(NPROCS 4)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 cube.msh -autoclean -partdual -metiskway ${NPROCS})

FILE(COPY ${BINARY_DIR}/fem/src/modules/ResultOutputSolve${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME ResultOutputSolve${SHLEXT} ResultOutputSolve1${SHLEXT})

RUN_ELMER_TEST()
