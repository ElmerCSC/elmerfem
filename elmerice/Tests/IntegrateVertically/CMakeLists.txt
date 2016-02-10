INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(integratevertically.sif integratevertically.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO mesh.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(IntegrateVertically)
ADD_ELMERICE_LABEL(IntegrateVertically elmerice-fast)
