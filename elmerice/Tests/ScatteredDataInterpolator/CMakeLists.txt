INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
FILE(COPY ELMERSOLVER_STARTINFO Rand200.txt FlightLines.txt mesh2D.grd README DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(ScatteredDataInterpolator)
ADD_ELMERICE_LABEL(ScatteredDataInterpolator netcdf)
