INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(ModelPDEsplitted ModelPDEsplitted ModelPDEsplitted.f90)

file(COPY square.grd ELMERSOLVER_STARTINFO square ModelPDEsplitted.f90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ModelPDEsplitted)
ADD_ELMER_LABEL(ModelPDEsplitted quick)
