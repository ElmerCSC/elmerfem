INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(BlockPoisson3 BlockPoissonAssembly BlockPoissonAssembly.f90)

file(COPY square.grd ELMERSOLVER_STARTINFO BlockPoissonAssembly.f90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(BlockPoisson3 LABELS block)
