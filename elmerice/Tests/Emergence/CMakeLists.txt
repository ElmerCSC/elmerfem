INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(emergencevelocity_restart.sif emergencevelocity_restart.sif COPYONLY)

ADD_ELMERICETEST_MODULE(Emergence buelerprofile buelerprofile.f90)
FILE(COPY ELMERSOLVER_STARTINFO buelerprofile.f90 footprint_bueler_f DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(Emergence)
ADD_ELMERICE_LABEL(Emergence elmerice-fast)
