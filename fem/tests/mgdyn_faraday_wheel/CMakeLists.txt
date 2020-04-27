INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(with_projection.sif with_projection.sif COPYONLY)

file(COPY Mesh ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_faraday_wheel LABELS mgdyn mortar)
