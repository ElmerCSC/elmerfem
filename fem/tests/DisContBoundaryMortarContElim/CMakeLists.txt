INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( discont.sif discont.sif COPYONLY)

file(COPY square.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(DisContBoundaryMortarContElim NPROCS 1 2 4 8 LABELS eliminate)
ADD_ELMER_LABEL(DisContBoundaryMortarContElim quick)
