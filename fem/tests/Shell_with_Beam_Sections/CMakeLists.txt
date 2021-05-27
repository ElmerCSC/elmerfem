INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( with_beams.sif with_beams.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO strip_with_beams DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Shell_with_Beam_Sections LABELS quick shell)
