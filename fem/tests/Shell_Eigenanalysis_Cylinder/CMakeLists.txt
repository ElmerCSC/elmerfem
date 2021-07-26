INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( eigenanalysis.sif eigenanalysis.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO fan_blade DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Shell_Eigenanalysis_Cylinder LABELS shell)
