IF(HAVE_ZOLTAN AND HAVE_MUMPS)
INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( eigenanalysis.sif eigenanalysis.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO cylinder DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Shell_with_Solid_EigenanalysisPar NPROCS 3 LABELS shell block eigen mumps zoltan)
ENDIF()
