INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( diffuser_sa.sif diffuser_sa.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO diffuser.grd diffuser DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(diffuser_sa)
ADD_ELMER_LABEL(diffuser_sa benchmark)
