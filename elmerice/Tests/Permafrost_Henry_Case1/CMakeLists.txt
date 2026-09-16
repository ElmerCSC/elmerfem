INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)

CONFIGURE_FILE(case.sif case.sif COPYONLY)
FILE(COPY ELMERSOLVER_STARTINFO material_case1.db solvent.db solute.db mesh
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMERICE_TEST(Permafrost_Henry_Case1)
ADD_ELMERICE_LABEL(Permafrost_Henry_Case1 permafrost)
