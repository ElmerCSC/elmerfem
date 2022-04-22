INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(ismip_weertman2.sif ismip_weertman2.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO rectangle.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(Friction_Weertman2)
ADD_ELMERICE_LABEL(Friction_Weertman2 elmerice-fast)
