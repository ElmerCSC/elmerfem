INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

FILE(COPY ${TEST_SOURCE}/../../../../install/share/elmersolver/lib/FreeSurfaceSolver.so DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME FreeSurfaceSolver.so MyFreeSurfaceSolver)

RUN_ELMERICE_TEST()
