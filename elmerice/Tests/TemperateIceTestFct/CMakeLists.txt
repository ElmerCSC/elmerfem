INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(temperaticetest_fct.sif temperaticetest_fct.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO column DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(TemperateIceTestFct)
ADD_ELMERICE_LABEL(TemperateIceTestFct elmerice-fast)
