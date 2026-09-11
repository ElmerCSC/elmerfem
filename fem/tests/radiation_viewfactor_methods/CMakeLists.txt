INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(numrays.sif numrays.sif COPYONLY)
CONFIGURE_FILE(cfrays.sif  cfrays.sif  COPYONLY)
CONFIGURE_FILE(cfclip.sif  cfclip.sif  COPYONLY)
CONFIGURE_FILE(numclip.sif numclip.sif COPYONLY)
CONFIGURE_FILE(cullon.sif  cullon.sif  COPYONLY)
CONFIGURE_FILE(culloff.sif culloff.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO box_in_box.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(radiation_viewfactor_methods LABELS radiation radiosity)
