INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( 1d.sif 1d.sif COPYONLY)
ADD_ELMERTEST_MODULE(1dtests_MeshLevels Poisson Poisson.f90)
ADD_ELMERTEST_MODULE(1dtests_MeshLevels Source Source.f90)

file(COPY Poisson.f90 ELMERSOLVER_STARTINFO Source.f90 1d.grd  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(1dtests_MeshLevels LABELS quick transient)
