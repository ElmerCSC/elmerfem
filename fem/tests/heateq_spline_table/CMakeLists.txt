INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( TempDist.sif TempDist.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO cond.dat angle.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(heateq_spline_table)
ADD_ELMER_LABEL(heateq_spline_table quick)
