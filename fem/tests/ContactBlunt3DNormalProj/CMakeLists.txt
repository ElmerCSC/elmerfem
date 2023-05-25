INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case.sif case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO base.grd blunt.msh extrude.eg DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ContactBlunt3DNormalProj LABELS contact mortar elasticsolve n-t)
