IF(HAVE_MUMPS)
INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(TEAM30a_3ph_w_bf.sif TEAM30a_3ph_w_bf.sif COPYONLY)
file(COPY ELMERSOLVER_STARTINFO TEAM30a_3ph DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(TEAM30a_3ph_w_bf LABELS quick harmonic)
ENDIF()
