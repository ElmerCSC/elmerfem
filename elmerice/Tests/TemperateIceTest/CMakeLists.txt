INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(temperaticetest.sif temperaticetest.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO column DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(TemperateIceTest)
ADD_ELMERICE_LABEL(TemperateIceTest elmerice-fast)
