INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( selfcapacitance.sif selfcapacitance.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO cubes.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(elstat_infty LABELS lumping)
