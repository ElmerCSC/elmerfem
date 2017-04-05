INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( reynolds_cylplate.sif reynolds_cylplate.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO cylplate.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(reynolds2)
ADD_ELMER_LABEL(reynolds2 quick)
