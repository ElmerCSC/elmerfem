INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(McKenzie.sif McKenzie.sif COPYONLY)

FILE(COPY McKenzie.db ELMERSOLVER_STARTINFO frozenwall DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(Permafrost_Frozenwall)
ADD_ELMERICE_LABEL(Permafrost_Frozenwall permafrost)
