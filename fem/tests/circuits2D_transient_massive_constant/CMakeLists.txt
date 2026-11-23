INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(sif/case.sif sif/case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY sif/1381-circuits-components.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
file(COPY 1381/mesh.boundary 1381/mesh.header 1381/mesh.nodes 1381/mesh.elements 1381/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1381")
file(COPY 1381/dat/1381.dat.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1381/dat")

IF(HAVE_UMFPACK OR HAVE_MUMPS)
  ADD_ELMER_TEST(circuits2D_transient_massive_constant NPROCS 1 4
    LABELS 2D circuits transient mgdyn quick)
ENDIF()
