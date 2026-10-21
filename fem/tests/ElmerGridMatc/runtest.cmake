include(test_macros)

if(WITH_MATC)
  execute_process(COMMAND ${ELMERGRID_BIN} 1 2 rect-matc -out rect)
else()
  execute_process(COMMAND ${ELMERGRID_BIN} 1 2 rect)
endif()

RUN_ELMER_TEST()
