INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( diffuser_sst.sif diffuser_sst.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO diffuser.grd diffuser DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(diffuser_sst)
