INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 PlanMesh.msh -autoclean -metis 2 0)

RUN_ELMERICE_TEST()
