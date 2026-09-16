INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
FILE(COPY ELMERSOLVER_STARTINFO lunardini_material.db lunardini_solvent.db box
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(Permafrost_Lunardini)
ADD_ELMERICE_LABEL(Permafrost_Lunardini permafrost)
