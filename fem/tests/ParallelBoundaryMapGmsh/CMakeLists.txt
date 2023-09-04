INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case_b.sif case_b.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO cube_b.grd case_a.msh_4np1 case_a.msh_4np2 case_a.msh_4np3 case_a.msh_4np4 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ParallelBoundaryMapGmsh LABELS restart NPROCS 1 4 8)
