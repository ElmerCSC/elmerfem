INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( TempDist.sif TempDist.sif COPYONLY)

file(COPY corona-mpi64.sh Mesh ELMERSOLVER_STARTINFO Mesh1.mif solver  mpd.hosts corona-mpi32.sh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(heateq-par)
ADD_ELMER_LABEL(heateq-par quick)
