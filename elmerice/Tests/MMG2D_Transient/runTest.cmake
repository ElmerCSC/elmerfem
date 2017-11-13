INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

FILE(COPY ${PROJECT_SOURCE_DIR}/..//Solvers/MeshAdaptation_2D/Script_Transient.sh  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

EXECUTE_PROCESS(COMMAND sh ./Script_Transient.sh)

RUN_ELMERICE_TEST()
