INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
CONFIGURE_FILE(case_singleiter.sif case_singleiter.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO Step.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Step_ns_vec_transient)
