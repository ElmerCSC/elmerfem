INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( esim1.sif esim1.sif COPYONLY)

file(COPY verkko.gid ELMERSOLVER_STARTINFO shape_disp.src  true_disp.src DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(stress)
ADD_ELMER_LABEL(stress quick)
