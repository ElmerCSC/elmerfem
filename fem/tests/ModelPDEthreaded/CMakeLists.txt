INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(ModelPDEthreaded ModelPDEthreaded ModelPDEthreaded.F90)

file(COPY square.grd ELMERSOLVER_STARTINFO ModelPDEthreaded.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ModelPDEthreaded)
ADD_ELMER_LABEL(ModelPDEthreaded threaded)
