IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(quad.sif  quad.sif  COPYONLY)
  CONFIGURE_FILE(tri.sif   tri.sif   COPYONLY)
  CONFIGURE_FILE(hex.sif   hex.sif   COPYONLY)
  CONFIGURE_FILE(prism.sif prism.sif COPYONLY)
  CONFIGURE_FILE(tet.sif   tet.sif   COPYONLY)

  file(COPY square.grd trimesh.grd cube.grd prismmesh.grd tet.msh ELMERSOLVER_STARTINFO
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(ElasticStabilized LABELS quick elasticity)
ENDIF()
