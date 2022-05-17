INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( openhemisphere.sif openhemisphere.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO hemisphere DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(Shell_OpenHemisphere_High_Order_Blending LABELS shell)
