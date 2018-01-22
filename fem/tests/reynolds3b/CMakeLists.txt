INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( reynolds_cylinder.sif reynolds_cylinder.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO cylinder.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(reynolds3b LABELS matc)
