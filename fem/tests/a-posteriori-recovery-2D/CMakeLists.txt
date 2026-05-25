INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( recovery_adapt.sif recovery_adapt.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO strip.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(a-posteriori-recovery-2D LABELS quick amr emwave)
