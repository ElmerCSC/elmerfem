INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(sif/case.sif sif/case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY sif/1381-circuits-components.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
file(COPY 1381/mesh.boundary 1381/mesh.header 1381/mesh.nodes 1381/mesh.elements 1381/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1381")
file(COPY 1381/dat/1381.dat.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1381/dat")

# The single-task run is kept so that the parallel result is compared against a
# serial one computed with the very same solver settings, not just against the
# stored reference norm.
ADD_ELMER_TEST(circuits2D_harmonic_stranded_parallel NPROCS 1 2
  LABELS 2D circuits harmonic mgdyn stranded)
