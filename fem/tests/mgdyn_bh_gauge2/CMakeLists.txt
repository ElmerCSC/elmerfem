INCLUDE(test_macros)
CONFIGURE_FILE(block.sif block.sif COPYONLY)
file(COPY ELMERSOLVER_STARTINFO mesh_tet HB DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
ADD_ELMER_TEST(mgdyn_bh_gauge2 LABELS whitney gauge block)
