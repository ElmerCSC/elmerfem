INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
FILE(COPY ELMERSOLVER_STARTINFO material.db solvent.db mesh
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(Permafrost_Frozenwall)
ADD_ELMERICE_LABEL(Permafrost_Frozenwall permafrost)
