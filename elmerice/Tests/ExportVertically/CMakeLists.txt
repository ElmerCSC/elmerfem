INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(exportvertically.sif exportvertically.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO mesh.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(ExportVertically)
ADD_ELMERICE_LABEL(ExportVertically elmerice-fast)
