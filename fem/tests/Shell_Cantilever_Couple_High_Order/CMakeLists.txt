INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( coupleload.sif coupleload.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO strip.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Shell_Cantilever_Couple_High_Order LABELS shell)
