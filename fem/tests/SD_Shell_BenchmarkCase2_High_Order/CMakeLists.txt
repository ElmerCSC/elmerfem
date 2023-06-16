INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( shell_cartesian.sif shell_cartesian.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO testmesh.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(SD_Shell_BenchmarkCase2_High_Order LABELS shell serendipity)
