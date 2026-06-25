INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(subshelfmeltlineartf.sif subshelfmeltlineartf.sif COPYONLY)

ADD_ELMERICETEST_MODULE(SubShelfMeltLinearTF InitOceanVarsSolver InitOceanVarsSolver.F90)

FILE(COPY ELMERSOLVER_STARTINFO README.txt mesh.grd InitOceanVarsSolver.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(SubShelfMeltLinearTF)
ADD_ELMERICE_LABEL(SubShelfMeltLinearTF elmerice-fast)
