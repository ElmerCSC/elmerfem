INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( flow_res_noslip.sif flow_res_noslip.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO hole.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(FlowResNoslip)
ADD_ELMER_LABEL(FlowResNoslip benchmark)
