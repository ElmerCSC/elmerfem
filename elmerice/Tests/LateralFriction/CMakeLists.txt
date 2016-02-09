INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(lateralFriction.sif lateralFriction.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO rectangle.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(LateralFriction)
ADD_ELMERICE_LABEL(LateralFriction elmerice-fast)
