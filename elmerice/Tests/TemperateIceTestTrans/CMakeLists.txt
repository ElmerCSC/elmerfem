INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(temperaticetest_transient.sif temperaticetest_transient.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO column DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(TemperateIceTestTrans)
ADD_ELMERICE_LABEL(TemperateIceTestTrans elmerice-fast)
