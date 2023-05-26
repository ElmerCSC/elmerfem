INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)

ADD_ELMERTEST_MODULE(AlternatingSource AlternatingSource AlternatingSource.F90)

file(COPY square.grd AlternatingSource.F90 ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(HeatControlExplicit LABELS quick control NPROCS 1 3)
