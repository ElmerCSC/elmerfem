INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(computeDevStress.sif computeDevStress.sif COPYONLY)

FILE(COPY ELMERSOLVER_STARTINFO rectangle.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(ComputeDevStress)
ADD_ELMERICE_LABEL(ComputeDevStress elmerice-fast)
