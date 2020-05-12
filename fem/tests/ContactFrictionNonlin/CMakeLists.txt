INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( nonlin.sif nonlin.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO blunt.msh blunt.geo DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ContactFrictionNonlin LABELS elasticsolve)
