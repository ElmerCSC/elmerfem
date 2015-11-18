INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

FILE(COPY ${BINARY_DIR}/fem/src/modules/FreeSurfaceSolver${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

FILE(RENAME FreeSurfaceSolver${SHLEXT}  MyFreeSurfaceSolver)

RUN_ELMERICE_TEST()
