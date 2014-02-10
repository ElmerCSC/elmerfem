
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 1d -decimals 20)
set(ENV{ELMER_LIB} "/home/silvonen/compile/ElmerGitSvn/code/fem/src")
set(ENV{ELMER_HOME} "/home/silvonen/compile/ElmerGitSvn/code/fem/src")
execute_process(COMMAND ${ELMERSOLVER_BIN} OUTPUT_FILE "test.log"
  ERROR_FILE "test.log")
execute_process(COMMAND ${FINDNORM_BIN} ${CMAKE_CURRENT_BINARY_DIR}/test.log
  OUTPUT_VARIABLE tulos)
message(STATUS ${tulos})