INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(Free_maxd.sif Free_maxd.sif COPYONLY)

ADD_ELMERTEST_MODULE(<testname> <modulename> <modulesource.f90>)

file(COPY Mesh ELMERSOLVER_STARTINFO Mesh.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(freesurf_maxd)

