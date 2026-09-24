INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case.sif case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO box_in_box.grd box_in_box_rad.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(radiation_box_in_box_twomesh NPROCS 1 2 4 LABELS radiation radiosity)
