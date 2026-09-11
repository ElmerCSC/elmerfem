INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( Step_sst-kw-wf_vec.sif Step_sst-kw-wf_vec.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO Step Step.grd  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Step_sst-kw-wf_vec)
ADD_ELMER_LABEL(Step_sst-kw-wf_vec benchmark)
