INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

SET(NPROCS 4)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 PlanMesh.msh -autoclean -metis ${NPROCS} 0)

#Calving3D depends on ElmerGrid - point to the just-built copy
SET(old_path "$ENV{PATH}")
IF(WIN32)
  SET(ENV{PATH} "${BINARY_DIR}/elmergrid/src;$ENV{PATH}")
ENDIF(WIN32)
IF(NOT(WIN32))
  SET(ENV{PATH} "${BINARY_DIR}/elmergrid/src:$ENV{PATH}")
ENDIF(NOT(WIN32))

RUN_ELMERICE_TEST(WITH_MPI)

#reset PATH
SET(ENV{PATH} "${old_path}")
