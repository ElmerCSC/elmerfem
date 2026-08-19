include(test_macros)

execute_process(COMMAND ${ELMERGRID_BIN} 1 2 rect.grd -partdual -metiskway ${MPIEXEC_NTASKS})

# Two routes through the redistribution of MUMPS's distributed solution, and both
# have to give the same answer. The gathered one goes first and its verdict is
# read here, because RUN_ELMER_TEST clears TEST.PASSED before running the sif
# named in ELMERSOLVER_STARTINFO, which would otherwise discard it.
EXECUTE_ELMER_SOLVER_MPI(gathered.sif)

SET(_gathered_result_file "TEST.PASSED_${MPIEXEC_NTASKS}")
IF(NOT EXISTS "${_gathered_result_file}")
  MESSAGE(FATAL_ERROR "gathered.sif wrote no result file")
ENDIF()
FILE(READ "${_gathered_result_file}" _gathered_result)
STRING(STRIP "${_gathered_result}" _gathered_result)
IF(NOT _gathered_result EQUAL "1")
  MESSAGE(FATAL_ERROR "gathered.sif failed its reference norm")
ENDIF()

RUN_ELMER_TEST()

# And assert which route the default run actually took. This is not belt and
# braces: a run that quietly fell back to the gathered redistribution looks
# exactly like one that used the direct mapping -- same norm, by design -- so
# without this the test would keep passing if the mapping stopped being reached,
# which is the one thing it exists to protect. The message is emitted at info
# level 8, hence the output level in case.sif.
FILE(READ "test-stdout_${MPIEXEC_NTASKS}.log" _default_run_output)
STRING(FIND "${_default_run_output}" "Solution redistribution plan built" _plan_used)
IF(_plan_used EQUAL -1)
  MESSAGE(FATAL_ERROR
    "case.sif did not use the distributed-solution mapping; it fell back to the "
    "gathered route, so this test verified nothing about it")
ENDIF()
