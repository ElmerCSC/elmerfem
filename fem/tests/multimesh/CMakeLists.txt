INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( multimesh.sif multimesh.sif COPYONLY)

file(COPY global.grd ELMERSOLVER_STARTINFO reference.grd local.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(multimesh LABELS quick)
