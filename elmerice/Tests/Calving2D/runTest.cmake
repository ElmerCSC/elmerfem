INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 block_helheim.msh -autoclean)


FILE(COPY ${BINARY_DIR}/fem/src/modules/FreeSurfaceSolver${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME FreeSurfaceSolver${SHLEXT} FreeSurfaceSolver1${SHLEXT})

FILE(COPY ${BINARY_DIR}/fem/src/modules/MeshSolve${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME MeshSolve${SHLEXT} MeshSolve1${SHLEXT})

RUN_ELMERICE_TEST()
