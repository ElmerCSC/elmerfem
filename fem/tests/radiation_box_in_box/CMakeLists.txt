INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case.sif case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO box_in_box.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(radiation_box_in_box LABELS radiation radiosity)
