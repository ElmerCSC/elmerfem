INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(calving2d.sif calving2d.sif COPYONLY)

ADD_ELMERICETEST_MODULE(Calving2D bedrockfunction_flat ./PROG/bedrockfunction_flat.f90)

FILE(COPY ELMERSOLVER_STARTINFO block_helheim.msh PROG DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(Calving2D)
ADD_ELMER_LABEL(Calving2D slow)
