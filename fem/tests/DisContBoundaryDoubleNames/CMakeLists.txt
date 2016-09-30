INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( discont.sif discont.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO circle_in_square DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(DisContBoundaryDoubleNames)
ADD_ELMER_LABEL(DisContBoundaryDoubleNames quick)
