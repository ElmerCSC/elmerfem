INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(ismip_coulomb_pw.sif ismip_coulomb_pw.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO rectangle.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(Friction_Coulomb_Pw)

