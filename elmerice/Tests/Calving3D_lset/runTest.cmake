INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

SET(NPROCS 2)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 PlanMesh.msh -autoclean -metis ${NPROCS} 0)

FILE(COPY ${BINARY_DIR}/fem/src/modules/FreeSurfaceSolver${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME FreeSurfaceSolver${SHLEXT} FreeSurfaceSolver1${SHLEXT})

FILE(COPY ${BINARY_DIR}/fem/src/modules/MeshSolve${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME MeshSolve${SHLEXT} MeshSolve1${SHLEXT})

FILE(COPY ${BINARY_DIR}/fem/src/modules/MeshSolve${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME MeshSolve${SHLEXT} MeshSolve2${SHLEXT})

FILE(COPY ${BINARY_DIR}/fem/src/modules/ResultOutputSolve${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME ResultOutputSolve${SHLEXT} ResultOutputSolve1${SHLEXT})

FILE(COPY ${BINARY_DIR}/fem/src/modules/ResultOutputSolve${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME ResultOutputSolve${SHLEXT} ResultOutputSolve2${SHLEXT})

FILE(COPY ${BINARY_DIR}/fem/src/modules/ResultOutputSolve${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME ResultOutputSolve${SHLEXT} ResultOutputSolve3${SHLEXT})

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
