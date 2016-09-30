INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( diffuser_v2f.sif diffuser_v2f.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO diffuser.grd diffuser DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(diffuser_v2f)
ADD_ELMER_LABEL(diffuser_v2f benchmark)
